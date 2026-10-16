"""Private importable entrypoint for interactive descriptor pools; trusted local transport only."""
from pathlib import Path
import pickle
import signal
import sys

from .runtime import native_threads, parallel_map


def _terminate(signum, frame):
    raise SystemExit(128 + signum)


def main():
    # Let the Pool context terminate/join workers when the parent cancels the broker.
    signal.signal(signal.SIGTERM, _terminate)
    root = Path(sys.argv[1])
    with root.joinpath('input').open('rb') as stream:
        function, work, count = pickle.load(stream)
    with native_threads():
        result = parallel_map(function, work, count)
    with root.joinpath('output').open('wb') as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
