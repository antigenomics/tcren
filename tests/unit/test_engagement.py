"""Small probability and held-out-data checks, independent of research inputs."""

import numpy as np
import polars as pl
from scipy.special import expit, logit


def test_empirical_tail_keeps_zero_and_amino_acid_capacity():
    from tcren.contacts.engagement import empirical_tails

    ref = pl.DataFrame(dict(aa=["A"] * 4 + ["R"] * 4,
                            contacted_atoms=[0, 1, 1, 2, 0, 2, 4, 6], weight=[1.] * 8))
    query = pl.DataFrame(dict(aa=["A", "A", "R", "A"], contacted_atoms=[0, 2, 2, 9]))
    bits = empirical_tails(ref, query)
    np.testing.assert_allclose(bits, -np.log2([1, 2/5, 4/5, 1/5]))
    assert bits[1] > bits[2]  # same count is rarer for alanine
    many = pl.DataFrame(dict(aa=["A"]*1000,contacted_atoms=[0,1]*500,weight=1/np.arange(1.,1001.)))
    np.testing.assert_array_equal(empirical_tails(many,query.filter(pl.col('contacted_atoms')==0)),[0.])


def test_offset_and_likelihood_only_use_training_labels():
    from tcren.contacts.engagement import offset_predict, likelihood_predict

    x = np.array([[0.], [0.], [1.], [1.], [0.], [1.]])
    y = np.array([0, 0, 1, 1, 0, 1])
    base = np.full(6, .5)
    train, test = np.arange(4), np.arange(4, 6)
    p, parts, params = offset_predict(x, y, base, train, test, .1)
    changed = y.copy(); changed[test] = 1 - changed[test]
    np.testing.assert_array_equal(p, offset_predict(x, changed, base, train, test, .1)[0])
    np.testing.assert_allclose(p, expit(logit(base[test]) + params["intercept"] + parts.sum(1)))
    assert p[0] < .5 < p[1]
    ranked = offset_predict(x,y,base,train,test,.1,loss="ap")[0]
    np.testing.assert_array_equal(ranked,offset_predict(x,changed,base,train,test,.1,loss="ap")[0])
    assert ranked[0] < .5 < ranked[1]
    for kind in ["states", "counts"]:
        n = np.ones_like(x)
        a = likelihood_predict(x, n, y, train, test, 4., kind)[0]
        b = likelihood_predict(x, n, changed, train, test, 4., kind)[0]
        np.testing.assert_array_equal(a, b)
        assert a[0] < .5 < a[1]


def test_transfer_keeps_n_terminal_anchors_and_c_terminal_alignment():
    from tcren.contacts.engagement import align_transfer

    names = ["tcr:TRA:1", "tcr:TRA:4", "tcr:TRA:5", "tcr:TRB:1", "tcr:TRB:4", "peptide:TRA:1"]
    source = ["tcr:TRA:1", "tcr:TRA:4", "tcr:TRA:5", "tcr:TRA:6",
              "tcr:TRB:1", "tcr:TRB:4", "tcr:TRB:5", "peptide:TRA:1"]
    x = np.arange(8.)[None, :]
    aligned, inventory = align_transfer({"joint": (np.empty((1,6)), names)},
                                        {"joint": (x, source)}, x+10)
    np.testing.assert_array_equal(aligned["joint"][0], [[0,2,3,4,6,7]])
    np.testing.assert_array_equal(inventory, [[10,12,13,14,16,17]])


def test_smooth_ap_gradient_and_order():
    from tcren.contacts.engagement import smooth_ap

    y = np.array([0., 1., 0., 1.])
    z = np.array([-.3, .1, .2, .8])
    ap, gradient = smooth_ap(z, y, .1)
    numerical = []
    for j in range(len(z)):
        delta = np.zeros(len(z)); delta[j] = 1e-6
        numerical.append((smooth_ap(z+delta,y,.1)[0]-smooth_ap(z-delta,y,.1)[0])/2e-6)
    np.testing.assert_allclose(gradient, numerical, atol=1e-8)
    assert 0 < ap < 1 and abs(gradient.sum()) < 1e-12
    assert smooth_ap(np.array([-10.,10.,-9.,11.]),y,.1)[0] > .999


def test_stable_selection_rejects_unstable_inner_gain_and_keeps_consistent_gain():
    from tcren.contacts.engagement import select_candidate

    search=[dict(E=.8,weight=0.,priority_AP=.5,worst_E=0.),
            dict(E=.9,weight=.5,priority_AP=.7,worst_E=-.02),
            dict(E=.85,weight=.2,priority_AP=.6,worst_E=.01)]
    assert select_candidate(search,'E',stable=True)['weight']==.2
    assert select_candidate(search,'E',stable=False)['weight']==.5


def test_training_quantile_transport_is_monotone_and_ignores_score_location():
    from tcren.contacts.engagement import transport_scores

    scores=np.array([.1,.2,.3,.4])
    baseline=np.array([.5,.6,.7,.8])
    query=np.array([.05,.15,.25,.45])
    mapped,params=transport_scores(scores,baseline,query)
    shifted,_=transport_scores(scores+.3,baseline,query+.3)
    np.testing.assert_allclose(mapped,shifted,atol=1e-15)
    assert np.all(np.diff(mapped)>=0) and mapped.min()>=.5 and mapped.max()<=.8
    np.testing.assert_array_equal(mapped,np.interp(query,params['score_knots'],params['baseline_knots']))


def test_shared_offset_dispatch_uses_training_labels_only():
    from tcren.contacts.engagement import prediction

    x=np.array([[0.],[0.],[1.],[1.],[0.],[1.]])
    matrices={'shared_core':(x,['contact_amount']),
              'baseline_missing':(np.zeros((6,1)),['baseline_missing'])}
    y=np.array([0,0,1,1,0,1]);train=np.arange(4);test=np.arange(4,6)
    a,_,params=prediction(('shared_core','shared_offset',.1),matrices,np.ones_like(x),y,np.full(6,.5),train,test)
    changed=y.copy();changed[test]=1-changed[test]
    b=prediction(('shared_core','shared_offset',.1),matrices,np.ones_like(x),changed,np.full(6,.5),train,test)[0]
    assert 'coef' in params and a[0]<.5<a[1]
    np.testing.assert_array_equal(a,b)


def test_bagged_prediction_averages_member_probabilities_without_test_labels():
    from tcren.contacts.engagement import prediction, bagged_prediction

    x=np.array([[0.],[1.],[0.],[1.],[0.],[1.],[0.],[1.]])
    matrices={'joint':(x,['amount']),'baseline_missing':(np.zeros((8,1)),['baseline_missing'])}
    y=np.array([0,1,0,1,0,0,0,1]);base=np.full(8,.5)
    trains=[np.arange(4),np.arange(2,6)];test=np.arange(6,8)
    spec=('joint','offset',.1)
    expected=np.mean([prediction(spec,matrices,np.ones_like(x),y,base,t,test)[0] for t in trains],axis=0)
    p,parts,params=bagged_prediction(spec,matrices,np.ones_like(x),y,base,trains,test,transport=False)
    np.testing.assert_allclose(p,expected,atol=1e-15)
    changed=y.copy();changed[test]=1-changed[test]
    q=bagged_prediction(spec,matrices,np.ones_like(x),changed,base,trains,test,transport=False)[0]
    np.testing.assert_array_equal(p,q)
    assert len(parts)==len(params['members'])==2
