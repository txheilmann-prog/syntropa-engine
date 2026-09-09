# Contributing to Syntropa

Thanks for your interest in improving Syntropa.

## Development setup

    python -m venv .venv
    . .venv/bin/activate            # Windows: .venv\Scripts\activate
    pip install -e ".[test]"
    pytest

## Reporting issues

Open an issue with a minimal reproducer: the component manifests or netlist, the command you ran, what you
expected, and what happened. Reports where a physics gate gives the wrong verdict (a composition that should be
rejected but is accepted, or the reverse) are especially valuable.

## Pull requests

- Keep each pull request focused on one concern.
- Add or update tests; `pytest` must pass.
- Model numbers must come from an authoritative, cited source, never fabricated; label nominal values as nominal.
- Run the formatter and linter (`ruff format`, `ruff check`) before submitting.

## Support

For questions, open a GitHub Discussion or an issue.
