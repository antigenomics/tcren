"""Watch a running producer and retain confidence before upstream removes its raw results."""

import argparse
import json
import pickle
from pathlib import Path
import time

from tcren.validation.tcrmodel2 import export_generated_confidence, retain_generated_inputs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=5)
    parser.add_argument("--retain-inputs", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    retained = {
        p.stem
        for p in args.output.glob("model_*.json")
        if json.loads(p.read_text()).get("token_order_verified_by_ca_plddt")
    }
    deadline = time.monotonic() + 8 * 3600
    last_status = None
    inputs_retained = (args.output/'producer-inputs.json').is_file()
    while True:
        if args.retain_inputs and not inputs_retained:
            try:
                retain_generated_inputs(args.input/'features.pkl',
                    args.input.parent/args.input.parent.name/'msa_features.pkl',args.output)
                inputs_retained = True
            except (ValueError, EOFError, FileNotFoundError, pickle.UnpicklingError):
                pass
        for path in sorted(args.input.glob("result_model_*.pkl")):
            model = path.stem.removeprefix("result_")
            coordinates = args.input / ("unrelaxed_" + model + ".pdb")
            output = args.output / (model + ".json")
            if model in retained or not coordinates.exists():
                continue
            # The producer writes the result pickle before writing its PDB. Both present means
            # confidence is complete; an incomplete PDB is retried on the next observation.
            try:
                export_generated_confidence(path, coordinates, output)
            except (ValueError, EOFError, FileNotFoundError):
                continue
            retained.add(model)
            print(f"Retained {model} confidence", flush=True)
        count = len(retained)
        state = dict(
            state="complete" if count == args.expected and (inputs_retained or not args.retain_inputs) else "waiting",
            files=count,
            expected=args.expected,
            source=str(args.input),
            inputs_retained=inputs_retained,
        )
        if state != last_status:
            (args.output / "retention-status.txt").write_text(json.dumps(state, indent=2))
            last_status = state
        if count == args.expected and (inputs_retained or not args.retain_inputs):
            break
        if time.monotonic() > deadline:
            raise TimeoutError("Raw confidence retention did not complete within 8 hours")
        time.sleep(0.05)
