"""Private-experiment building blocks for residue engagement ranking.

The fixed baseline is ipTM + recognition Betti topology. Contact additions are
selected using grouped inner predictions, never evaluation labels. All quantities
are ranking scores; neither mock tails nor logistic outputs imply affinity.
"""

import json
from pathlib import Path
import time
import warnings

import numpy as np
import polars as pl
from scipy.optimize import minimize
from scipy.special import expit, logit
import typer

from . import profile

MIX = tuple(np.round(np.linspace(0, 1, 21), 2))


def empirical_tails(reference, query):
    """Weighted type-specific inclusive survival, plus one tail pseudocount.

Weights have mean one among reference residues of each type. Consequently zero
has survival one, unseen large counts remain finite, and no atom-count division
or conditional-on-contact truncation is performed.
    """
    result = np.empty(query.height)
    for aa in query["aa"].unique():
        ref = reference.filter(pl.col("aa") == aa)
        if ref.is_empty():
            raise ValueError(f"No reference observations for amino acid {aa}")
        k = ref["contacted_atoms"].to_numpy()
        w = ref["weight"].to_numpy().astype(float)
        if not np.isfinite(w).all() or (w <= 0).any() or (k < 0).any():
            raise ValueError("Reference weights/counts must be finite and positive/nonnegative")
        w = w / w.mean()
        order = np.argsort(k, kind="stable")
        cumulative = np.r_[np.cumsum(w[order][::-1])[::-1], 0.]
        take = (query["aa"] == aa).to_numpy()
        counts = query["contacted_atoms"].to_numpy()[take]
        if not np.isfinite(counts).all() or (counts < 0).any():
            raise ValueError("Query counts must be finite and nonnegative")
        indices = np.searchsorted(k[order], counts, side="left")
        # Use the identical accumulated total in numerator/denominator: S(0)=1
        # must be exact, or standardization amplifies roundoff into fake signal.
        result[take] = -np.log2((1 + cumulative[indices]) / (1 + cumulative[0]))
    return result


def smooth_ap(scores, y, temperature=.1):
    """Full-cohort smooth AP and its analytic score gradient; ties use sigmoid.

The temporary matrix has n_positive by n_cohort elements, restricted here to
one small epitope cohort. Self-comparisons contribute one, with zero derivative.
    """
    positive = np.flatnonzero(y == 1)
    if not len(positive) or len(positive) == len(y) or temperature <= 0:
        raise ValueError("Smooth AP needs both classes and positive temperature")
    pair = expit((scores[None,:]-scores[positive,None])/temperature)
    pair[np.arange(len(positive)),positive] = 0.
    numerator = 1+pair @ y
    denominator = 1+pair.sum(1)
    influence = ((y[None,:]*denominator[:,None]-numerator[:,None])
                 /denominator[:,None]**2/len(positive))
    influence *= pair*(1-pair)/temperature
    gradient = influence.sum(0)
    gradient[positive] -= influence.sum(1)
    return float(np.mean(numerator/denominator)), gradient


def offset_predict(x, y, baseline, train, test, c, *, loss="logistic"):
    """Ridge conditional log-odds correction to a frozen baseline ranking."""
    from sklearn.preprocessing import StandardScaler

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nan_to_num(np.nanmean(x[train], axis=0), nan=0.)

    def design(indices):
        return np.c_[np.where(np.isfinite(x[indices]), x[indices], mean),
                     ~np.isfinite(x[indices])].astype(float)

    a, b = design(train), design(test)
    scaler = StandardScaler().fit(a)
    center, scale = scaler.mean_, scaler.scale_
    a, b = scaler.transform(a), scaler.transform(b)
    offset = logit(np.clip(baseline, 1e-6, 1-1e-6))

    def objective(theta):
        z = offset[train] + theta[0] + a @ theta[1:]
        if loss == "ap":
            ap, gradient = smooth_ap(z, y[train])
            residual = -len(train)*gradient
            cost = len(train)*(1-ap)
        else:
            residual = expit(z) - y[train]
            cost = np.logaddexp(0,z).sum()-y[train] @ z
        return (cost + .5/c * (theta[1:] @ theta[1:]),
                np.r_[residual.sum(), a.T @ residual + theta[1:]/c])

    fit = minimize(objective, np.zeros(a.shape[1]+1), jac=True, method="L-BFGS-B",
                   options={"maxiter": 1000, "gtol": 1e-6, "ftol": 1e-12})
    if not fit.success:
        raise RuntimeError(f"Offset logistic did not converge: {fit.message}")
    parts = b * fit.x[1:]
    params = dict(imputation=mean.tolist(), mean=center.tolist(), scale=scale.tolist(),
                  coef=fit.x[1:].tolist(), intercept=float(fit.x[0]),loss=loss)
    return expit(offset[test]+fit.x[0]+parts.sum(1)), parts, params


def likelihood_predict(x, available, y, train, test, strength, kind):
    """Signed smoothed per-position likelihood ratios, including gaps as a state.

Three engagement states (zero/partial/full) plus gaps, or exact integer counts
0..20 plus gaps. A symmetric Dirichlet prior supplies total mass `strength`.
The summed log ratio is a composite likelihood, not independent evidence.
    """
    missing = ~np.isfinite(x)
    if kind == "states":
        state = np.where(x == 0, 0, np.where(x == available, 2, 1))
        levels = 4
    elif kind == "counts":
        if np.any(np.isfinite(x) & ((x < 0) | (x > 20) | (x != np.floor(x)))):
            raise ValueError("Heavy-atom counts must be integers in 0..20")
        state = np.nan_to_num(x).astype(int)
        levels = 22
    else:
        raise ValueError("Unknown likelihood representation")
    state[missing] = levels-1
    table = []
    for label in [0, 1]:
        observed = state[train[y[train] == label]]
        counts = (observed[:, :, None] == np.arange(levels)).sum(0)
        table.append((counts + strength/levels) / (len(observed)+strength))
    llr = np.log(table[1]/table[0])
    parts = llr[np.arange(x.shape[1])[None, :], state[test]]
    return expit(parts.sum(1)/4.), parts/4., dict(log_ratio=llr.tolist(), temperature=4.)


def feature_matrices(models, detail, original, reference):
    """Retain existing raw profiles and add empirical type tails and state laws."""
    matrices = {k: original[k] for k in ["commitment", "joint", "alpha", "peptide", "states"]}
    prof = detail.filter((pl.col("side") == "peptide") | (pl.col("region") == "CDR3"))
    prof = prof.with_columns(pl.concat_str("side", "chain", "aligned_position", separator=":").alias("feature"))
    prof = prof.with_columns(pl.Series("tail_bits", empirical_tails(reference, prof)))
    wide = prof.select("id", "feature", "tail_bits").pivot(on="feature", index="id", values="tail_bits")
    names = sorted(c for c in wide.columns if c != "id")
    matrices["tail_profile"] = profile.ordered_matrix(models, wide, names), names
    sums = (prof.group_by("id", "side", "chain").agg(pl.col("tail_bits").sum())
            .with_columns(pl.concat_str("side", "chain", separator=":").alias("feature"))
            .pivot(on="feature", index="id", values="tail_bits"))
    names = sorted(c for c in sums.columns if c not in ["id", "side", "chain"])
    matrices["tail_sum"] = profile.ordered_matrix(models, sums, names), names
    for family in ["commitment", "tail_sum"]:
        x, names = matrices[family]
        contrast = models["iptm"].to_numpy()-models["prepared"].to_numpy()
        matrices[family+"_confidence"] = np.c_[x,contrast], names+["confidence_contrast"]
    x,names=matrices["commitment"]
    chosen=[j for j,name in enumerate(names) if name.split(":")[-1] in
            ["contacted_atoms","fully_contacted","involved_residues"]]
    tail,tail_names=matrices["tail_sum"]
    contrast=models["iptm"].to_numpy()-models["prepared"].to_numpy()
    matrices["shared_core"]=(np.c_[x[:,chosen],tail,contrast],
        [names[j] for j in chosen]+["tail:"+n for n in tail_names]+["confidence_contrast"])
    matrices["baseline_missing"] = (~np.isfinite(models["prepared"].to_numpy().astype(float)))[:,None], ["baseline_missing"]
    wide = prof.select("id", "feature", "available_atoms").pivot(on="feature", index="id", values="available_atoms")
    count_names = matrices["joint"][1]
    available = profile.ordered_matrix(models, wide, count_names)
    return matrices, available


def align_transfer(training, transfer, available):
    """Align C-terminal CDR3 padding to training columns; ignore novel insertions."""
    shifts = {}
    for chain in ["TRA", "TRB"]:
        prefix = f"tcr:{chain}:"
        def maximum(names):
            return max(int(n.split(":")[2]) for n in names if n.startswith(prefix))
        shifts[chain] = maximum(transfer["joint"][1])-maximum(training["joint"][1])

    def rename(name):
        parts = name.split(":")
        if parts[0] == "tcr" and len(parts) >= 3 and parts[2].isdigit() and int(parts[2]) > 3:
            position = int(parts[2])-shifts[parts[1]]
            if position <= 3:
                return None
            parts[2] = str(position)
        return ":".join(parts)

    aligned = {}
    for family, (_, names) in training.items():
        x, source_names = transfer[family]
        lookup = {rename(n):j for j,n in enumerate(source_names) if rename(n) is not None}
        dest = np.full((len(x), len(names)), np.nan)
        for j,name in enumerate(names):
            if name in lookup:
                dest[:,j] = x[:,lookup[name]]
        aligned[family] = dest, names
    lookup = {rename(n):j for j,n in enumerate(transfer["joint"][1]) if rename(n) is not None}
    aligned_available = np.full_like(aligned["joint"][0], np.nan)
    for j,name in enumerate(training["joint"][1]):
        if name in lookup:
            aligned_available[:,j] = available[:,lookup[name]]
    return aligned, aligned_available


def metric_vector(models, values):
    from sklearn.metrics import average_precision_score

    y, epi = models["y"].to_numpy(), models["cohort"].to_numpy()
    return np.array([average_precision_score(y[epi == e], values[epi == e]) for e in sorted(set(epi))])


def select_candidate(search, epitope, *, stable=False):
    """Choose inner AP; optionally require no loss in any inner fold."""
    eligible = [r for r in search if not stable or r["worst_"+epitope] >= 0]
    return max(eligible, key=lambda r:(r[epitope], -r["weight"], r["priority_AP"]))


def candidates(matrices, ranking=False, shared=False):
    result = ([(name, mode, c) for name in matrices if not name.endswith("_confidence") and name not in ["baseline_missing","shared_core"]
               for mode in ["logistic", "offset"] for c in profile.CS]
              + [("joint", mode, a) for mode in ["state_llr", "count_llr"] for a in [4., 16., 64.]])
    if ranking:
        result += [(name,"ap",c) for name in ["joint","commitment_confidence","tail_sum_confidence"] for c in profile.CS]
    if shared:
        result += [("shared_core",mode,c) for mode in ["shared_logistic","shared_offset"] for c in profile.CS]
    return result


def transport_scores(training_scores, training_baseline, scores):
    """Map learned scores to baseline units using training empirical quantiles.

Midrank percentiles handle ties. No query labels or query-cohort distribution
enters the mapping; bounded baseline quantiles give a bounded monotone map.
    """
    if not len(training_scores) or len(training_scores) != len(training_baseline):
        raise ValueError("Transport needs paired nonempty training scores/baselines")
    if not all(np.isfinite(a).all() for a in [training_scores,training_baseline,scores]):
        raise ValueError("Transport needs finite values")
    if np.any((training_baseline < 0) | (training_baseline > 1)):
        raise ValueError("Transport baseline must be in [0,1]")
    knots, counts = np.unique(training_scores, return_counts=True)
    quantiles = (np.cumsum(counts)-counts/2)/counts.sum()
    values = np.quantile(training_baseline,quantiles)
    return np.interp(scores,knots,values), dict(score_knots=knots.tolist(),baseline_knots=values.tolist())


def prediction(spec, matrices, available, y, base, train, test, *, transport=False):
    if transport:
        test = np.r_[train,test]
    name, mode, c = spec
    mode=mode.removeprefix("shared_")
    x = matrices[name][0]
    names = matrices[name][1]
    if mode in ["logistic", "offset", "ap"]:
        x = np.c_[x,matrices["baseline_missing"][0]]
        names = names+["baseline_missing"]
    if mode == "logistic":
        result = profile.fit_predict(x, y, train, test, c)
    elif mode in ["offset", "ap"]:
        result = offset_predict(x, y, base, train, test, c, loss="ap" if mode == "ap" else "logistic")
    else:
        result = likelihood_predict(x, available, y, train, test, c,
                                    "states" if mode == "state_llr" else "counts")
    result[2]["features"] = names+["missing:"+n for n in names] if mode in ["logistic","offset","ap"] else names
    if transport:
        score,parts,params = result
        mapped,mapping = transport_scores(score[:len(train)],base[train],score[len(train):])
        params["transport"] = mapping
        result = mapped,parts[len(train):],params
    return result


def bagged_prediction(spec, matrices, available, y, raw_base, training_sets, test, *, transport=True):
    """Average probabilities from the already-defined inner training partitions."""
    scores, parts, members = [], [], []
    for train in training_sets:
        if set(y[train]) != {0,1}:
            raise ValueError("Every averaging member needs both training classes")
        imputation = float(np.nanmean(raw_base[train]))
        base = np.where(np.isfinite(raw_base),raw_base,imputation)
        score, contribution, params = prediction(spec,matrices,available,y,base,train,test,transport=transport)
        params["baseline_imputation"] = imputation
        scores.append(score);parts.append(contribution);members.append(params)
    return np.mean(scores,axis=0),parts,dict(members=members)


def fold_job(task):
    """One coarse fold job; native kernels are limited by runtime.parallel_map."""
    models, matrices, available, fold, train, test, policy, transfer = task
    y = models["y"].to_numpy()
    epi = models["cohort"].to_numpy()
    # Missing baseline: use training-only mean, and keep the established mask.
    raw_base = models["prepared"].to_numpy().astype(float)
    base = np.where(np.isfinite(raw_base), raw_base, np.nanmean(raw_base[train]))
    dev = models[train]
    inner = profile.splits(dev, 3, profile.SEED+100+fold)
    epochs = sorted(set(epi))
    priority = [epochs.index(e) for e in profile.PRIORITY]
    baseline_ap = metric_vector(dev, base[train])
    baseline_inner_ap = np.array([metric_vector(dev[iv],base[train[iv]]) for _,iv in inner])
    choices = [("baseline", "baseline", 0.)] + candidates(matrices, ranking=policy in ["ranking","stable","transport","shared","bagged"],shared=policy in ["shared","bagged"])
    predictions = []
    search = []
    for index, spec in enumerate(choices):
        pred = base[train].copy()
        if index:
            for it, iv in inner:
                if spec[1].startswith("shared_"):
                    a,b=train[it],train[iv]
                    inner_base=np.where(np.isfinite(raw_base),raw_base,np.nanmean(raw_base[a]))
                    pred[iv]=prediction(spec,matrices,available,y,inner_base,a,b,transport=True)[0]
                    continue
                for e in epochs:
                    a = train[it[epi[train[it]] == e]]
                    b = train[iv[epi[train[iv]] == e]]
                    # Imputation of baseline is confined to this inner training set.
                    inner_base = np.where(np.isfinite(raw_base), raw_base, np.nanmean(raw_base[a]))
                    pred[iv[epi[train[iv]] == e]] = prediction(spec, matrices, available, y, inner_base, a, b,transport=policy in ["transport","shared","bagged"])[0]
        predictions.append(pred)
        for weight in ([0.] if index == 0 else MIX[1:]):
            mix = (1-weight)*base[train] + weight*pred
            ap = metric_vector(dev, mix)
            fold_delta = np.array([metric_vector(dev[iv],mix[iv]) for _,iv in inner])-baseline_inner_ap
            search.append(dict(candidate=index, family=spec[0], method=spec[1], regularization=spec[2],
                               weight=float(weight), priority_AP=float(ap[priority].mean()),
                               macro_AP=float(ap.mean()), worst_delta=float((ap-baseline_ap).min()),
                               **{"worst_"+e:float(v) for e,v in zip(epochs,fold_delta.min(0))},
                               **{e:float(v) for e,v in zip(epochs, ap)}))
    # Global priority selection is retained. Safeguard makes cohort-specific choices
    # using only inner predictions, with baseline as an explicit candidate.
    best_global = max(search, key=lambda r:(r["priority_AP"], r["macro_AP"], -r["weight"]))
    selected = []
    result = base[test].copy()
    global_result = base[test].copy()
    attributions = []
    for e in epochs:
        best = best_global if policy == "priority" else select_candidate(search,e,stable=policy == "stable")
        a, b = train[epi[train] == e], test[epi[test] == e]
        params = None
        priority_params = None
        for choice, output in [(best, result), (best_global, global_result)]:
            spec = choices[choice["candidate"]]
            if choice["candidate"]:
                fit_train=train if spec[1].startswith("shared_") else a
                if policy == "bagged":
                    training_sets=[train[it] if spec[1].startswith("shared_") else train[it][epi[train[it]] == e] for it,_ in inner]
                    p,part_members,fitted=bagged_prediction(spec,matrices,available,y,raw_base,training_sets,b)
                    member_params=fitted["members"]
                else:
                    p, parts, fitted = prediction(spec, matrices, available, y, base, fit_train, b,transport=policy in ["transport","shared"])
                    part_members,member_params=[parts],[fitted]
                output[epi[test] == e] = (1-choice["weight"])*base[b] + choice["weight"]*p
                if output is result:
                    params = fitted
                    for member,(parts,par) in enumerate(zip(part_members,member_params)):
                        contribution = (pl.DataFrame(parts,schema=par["features"])
                            .with_columns(pl.Series("id",models[b]["id"]))
                            .unpivot(index="id",variable_name="feature",value_name="logit_contribution")
                            .with_columns(pl.lit(fold).alias("fold"),pl.lit(e).alias("cohort"),
                                          pl.lit(member).alias("member"),
                                          pl.lit(float(choice["weight"])).alias("contact_weight")))
                        attributions.append(contribution)
                else:
                    priority_params = fitted
        selected.append(dict(fold=fold, cohort=e, **best, parameters=params,
                             priority_model=dict(**best_global,parameters=priority_params)))
    rows = models[test].select("id", "cohort", "y", "matched", "group").with_columns(
        pl.lit(fold).alias("fold"), pl.Series("baseline", base[test]),
        pl.Series("selected", result), pl.Series("priority_selected", global_result))
    transfer_rows = None
    transfer_log = []
    if transfer is not None:
        combined_models, combined_matrices, combined_available = transfer
        n = models.height
        target = combined_models[n:]
        both_y = combined_models["y"].to_numpy()
        both_raw_base = combined_models["prepared"].to_numpy().astype(float)
        both_base = both_raw_base.copy()
        both_base = np.where(np.isfinite(both_base), both_base, np.nanmean(raw_base[train]))
        forbidden_a, forbidden_b = set(target["cdr3_alpha"]), set(target["cdr3_beta"])
        allowed = np.array([(a not in forbidden_a and b not in forbidden_b)
                            for a,b in zip(models["cdr3_alpha"], models["cdr3_beta"])])
        safe_train = train[allowed[train]]
        output = both_base[n:].copy()
        for (e,), part in target.with_row_index("_index").group_by("cohort"):
            target_indices = part["_index"].to_numpy()
            b = n+target_indices
            known = e in epochs
            a = safe_train[epi[safe_train] == e] if known else safe_train
            if set(both_y[a]) != {0,1}:
                raise ValueError(f"Transfer training lacks two classes: {e}")
            best = next(s for s in selected if s["cohort"] == e) if known else best_global
            if best["candidate"]:
                if choices[best["candidate"]][1].startswith("shared_"):
                    a=safe_train
                if policy == "bagged":
                    training_sets=[]
                    for it,_ in inner:
                        member=train[it][allowed[train[it]]]
                        if known and not choices[best["candidate"]][1].startswith("shared_"):
                            member=member[epi[member] == e]
                        training_sets.append(member)
                    pred=bagged_prediction(choices[best["candidate"]],combined_matrices,combined_available,
                                           both_y,both_raw_base,training_sets,b)[0]
                else:
                    pred = prediction(choices[best["candidate"]], combined_matrices, combined_available,
                                      both_y, both_base, a, b,transport=policy in ["transport","shared"])[0]
                output[target_indices] = (1-best["weight"])*both_base[b]+best["weight"]*pred
            transfer_log.append(dict(fold=fold,cohort=e,known_epitope=known,train_n=len(a),
                                     excluded_shared_chain_n=int((~allowed[train]).sum()),
                                     family=best["family"],method=best["method"],weight=best["weight"]))
        transfer_rows = target.select("id", "cohort", "stratum", "y", "matched").with_columns(
            pl.lit(fold).alias("fold"), pl.Series("baseline",both_base[n:]), pl.Series("selected",output))
    return rows, selected, [dict(fold=fold, **s) for s in search], transfer_rows, transfer_log, attributions


def report_metrics(rows, out):
    from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

    records, curves = [], []
    for (e,), p in rows.group_by("cohort"):
        for population, q in [("all_labelled", p), ("observed_baselines", p.filter("matched"))]:
            y = q["y"].to_numpy()
            ref = average_precision_score(y, q["baseline"])
            for arm in ["baseline", "selected", "priority_selected"]:
                values = q[arm].to_numpy()
                if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
                    raise ValueError("Nonfinite/unbounded score")
                ap = average_precision_score(y, values)
                precision, recall, thresholds = precision_recall_curve(y, values)
                records.append(dict(cohort=e, population=population, arm=arm, n=q.height,
                                    binders=int(y.sum()), nonbinders=int((1-y).sum()), AP=float(ap),
                                    delta_AP=float(ap-ref), AUROC=float(roc_auc_score(y, values))))
                curves.extend(dict(cohort=e, population=population, arm=arm, precision=float(pr),
                                   recall=float(re), threshold=float(th))
                              for pr,re,th in zip(precision, recall, np.r_[thresholds, np.nan]))
    table = pl.DataFrame(records).sort("population", "cohort", "arm")
    table.write_csv(out/"metrics.tsv", separator="\t")
    pl.DataFrame(curves).write_csv(out/"precision_recall.tsv", separator="\t")
    return table


def command(
    contacts: Path = typer.Option(..., exists=True),
    atoms: Path = typer.Option(..., exists=True),
    baseline: Path = typer.Option(..., exists=True),
    mock_contacts: Path = typer.Option(..., exists=True),
    mock_atoms: Path = typer.Option(..., exists=True),
    mock_baseline: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    policy: str = typer.Option("priority"),
    threads: int = typer.Option(0, "--threads", "-t"),
):
    """Compare contact additions to a frozen ipTM+Betti baseline in grouped folds."""
    from sklearn.metrics import average_precision_score, roc_auc_score

    from tcren.runtime import parallel_map, resolve_threads

    if policy not in ["priority", "safeguard", "ranking", "stable", "transport", "shared", "bagged"]:
        raise ValueError("Unknown engagement policy")
    budget = resolve_threads(threads)
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    raw = output/"query"; raw.mkdir()
    mock = output/"reference"; mock.mkdir()
    all_models, all_original = profile.prepare(contacts, atoms, baseline, raw, include_unlabelled=True)
    mock_models, _ = profile.prepare(mock_contacts, mock_atoms, mock_baseline, mock)
    if not (mock_models["y"] == 0).all():
        raise ValueError("Mock reference requires only nonbinder/mock outcomes")
    for column in ["iptm","prepared"]:
        values=all_models[column].to_numpy().astype(float)
        if np.isinf(values).any() or np.any(np.isfinite(values) & ((values < 0) | (values > 1))):
            raise ValueError(f"{column} must be missing or bounded in [0,1]")
    # The external reference must not carry evaluation epitopes or either chain.
    forbidden_a = set(all_models["cdr3_alpha"])
    forbidden_b = set(all_models["cdr3_beta"])
    forbidden_e = set(all_models["cohort"])
    keep = mock_models.filter(~pl.col("cdr3_alpha").is_in(forbidden_a)
                             & ~pl.col("cdr3_beta").is_in(forbidden_b)
                             & ~pl.col("cohort").is_in(forbidden_e))
    if keep.is_empty():
        raise ValueError("No disjoint mock reference remains")
    ref = pl.read_parquet(mock/"residues.parquet").join(keep.select("id"), on="id", how="semi")
    ref = ref.with_columns(
        pl.when(pl.col("chain") == "TRA").then(pl.col("cdr3_alpha")).otherwise(pl.col("cdr3_beta")).alias("receptor"))
    ref = ref.with_columns(
        (1/pl.col("id").n_unique().over("cohort", "side", "chain", "region", "receptor")
         /pl.col("receptor").n_unique().over("cohort", "side", "chain", "region")).alias("weight"))
    # Normalize equal epitope/CDR weights to mean one before pseudo-count smoothing.
    ref = ref.with_columns((pl.col("weight")/pl.col("weight").mean()).alias("weight"))
    ref.write_parquet(output/"mock_reference.parquet")
    keep.write_parquet(output/"mock_models.parquet")
    models = all_models.filter(pl.col("dataset") == "tcrvdb").sort("id")
    models = models.with_columns(pl.Series("group", profile.connected_groups(models["cdr3_alpha"], models["cdr3_beta"])))
    # Reuse the same structure-to-profile implementation on TCR rows to preserve
    # the previous position padding and exact-chain fold membership.
    tcr_baseline = output/"tcr_baseline.parquet"
    pl.read_parquet(baseline).filter(pl.col("dataset") == "tcrvdb").write_parquet(tcr_baseline)
    tcr = output/"tcr"; tcr.mkdir()
    models, original = profile.prepare(contacts, atoms, tcr_baseline, tcr)
    detail = pl.read_parquet(tcr/"residues.parquet")
    matrices, available = feature_matrices(models, detail, original, ref)
    all_matrices, all_available = feature_matrices(all_models, pl.read_parquet(raw/"residues.parquet"), all_original, ref)
    aligned, aligned_available = align_transfer(matrices, all_matrices, all_available)
    lookup = {s:i for i,s in enumerate(all_models["id"])}
    indices = np.array([lookup[s] for s in models["id"]])
    for family in matrices:
        np.testing.assert_allclose(aligned[family][0][indices], matrices[family][0], atol=1e-12, equal_nan=True)
    np.testing.assert_allclose(aligned_available[indices], available, equal_nan=True)
    target_indices = np.array([i for i,s in enumerate(all_models["dataset"]) if s == "vdjdb"])
    order = np.r_[indices, target_indices]
    transfer = (all_models[order], {k:(v[0][order],v[1]) for k,v in aligned.items()}, aligned_available[order])
    folds = profile.splits(models, 5, profile.SEED)
    results = parallel_map(fold_job, [(models, matrices, available, i, a, b, policy, transfer)
                                     for i,(a,b) in enumerate(folds)], budget)
    rows = pl.concat([r[0] for r in results]).sort("id")
    if rows.height != models.height or rows["id"].n_unique() != models.height:
        raise ValueError("OOF row conservation failed")
    rows.write_parquet(output/"oof_predictions.parquet")
    table = report_metrics(rows, output)
    (output/"selected.json").write_text(json.dumps([s for r in results for s in r[1]], indent=2)+"\n")
    pl.DataFrame([s for r in results for s in r[2]]).write_csv(output/"inner_search.tsv", separator="\t")
    if any(r[5] for r in results):
        pl.concat([s for r in results for s in r[5]]).write_parquet(output/"feature_contributions.parquet")
    transfer_rows = pl.concat([r[3] for r in results])
    transfer_rows.write_parquet(output/"transfer_fold_predictions.parquet")
    pl.DataFrame([s for r in results for s in r[4]]).write_csv(output/"transfer_training.tsv",separator="\t")
    transfer_rows = transfer_rows.group_by("id","cohort","stratum","y","matched").agg(pl.col("baseline","selected").mean())
    transfer_rows.write_parquet(output/"transfer_predictions.parquet")
    reference_metrics = []
    for (stratum,e), part in transfer_rows.filter("matched").group_by("stratum","cohort"):
        y = part["y"].to_numpy()
        ref_ap = average_precision_score(y,part["baseline"])
        for arm in ["baseline","selected"]:
            ap = average_precision_score(y,part[arm])
            reference_metrics.append(dict(stratum=stratum,cohort=e,arm=arm,n=part.height,
                binders=int(y.sum()),nonbinders=int((1-y).sum()),AP=float(ap),delta_AP=float(ap-ref_ap),
                AUROC=float(roc_auc_score(y,part[arm]))))
    pl.DataFrame(reference_metrics).sort("stratum","cohort","arm").write_csv(output/"reference_metrics.tsv",separator="\t")
    (output/"PROVENANCE.json").write_text(json.dumps(dict(seed=profile.SEED, policy=policy,
        baseline="Fixed prior ipTM+Betti prepared score; no independent ipTM re-addition",
        models=models.height, groups=models["group"].n_unique(), mock_models=keep.height,
        mock_epitopes=keep["cohort"].n_unique(), workers=min(budget,5),
        matrix_shapes={k:list(v[0].shape) for k,v in matrices.items()}, seconds=time.monotonic()-start,
        inputs={str(p):profile.digest(p) for p in [contacts,atoms,baseline,mock_contacts,mock_atoms,mock_baseline]},
        interpretation="Exploratory internal OOF after prior panel inspection; no external validation or default promotion"), indent=2)+"\n")
    print(table.filter((pl.col("population") == "observed_baselines") & (pl.col("arm") != "priority_selected")), flush=True)
