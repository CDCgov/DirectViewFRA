# Contributing to DirectViewFRA

Thank you for contributing. Contributions that improve correctness, usability, documentation, testing, and reproducibility are welcome.

Before contributing, read the [`README`](README.md), [`LICENSE`](LICENSE), [`DISCLAIMER`](DISCLAIMER.md), and [`code of conduct`](code-of-conduct.md).

## Protect data and confidentiality

Only non-sensitive, publicly releasable information may be submitted. Do not include personally identifiable information, protected health information, non-public laboratory identifiers, credentials, internal infrastructure details, or restricted data in issues, pull requests, screenshots, logs, or example files.

Use synthetic or fully cleared example data when reproducing a problem.

## Report a problem

Open a GitHub issue and include:

- the operating system;
- R and package versions;
- whether hardcoded or command-line settings were used;
- the relevant settings, with sensitive paths and identifiers removed;
- the exact error or unexpected behavior; and
- a minimal, non-sensitive example when possible.

## Submit a pull request

A pull request should:

- explain the analytical or usability problem addressed;
- preserve the documented output schema unless a schema change is intentional and documented;
- add or update tests for changed behavior;
- update user documentation when settings, inputs, outputs, or reporting rules change;
- avoid duplicate calculated fields or undocumented fallback behavior;
- preserve Windows and Linux path compatibility; and
- confirm that no sensitive information is included.

For changes to normalization, NT50 calculation, reportability, or quality-control behavior, describe the scientific rationale and include representative test cases for well-behaved, non-bracketing, irregular, and failed-fit curves as applicable.

## Public domain

This project is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). All contributions will be released under the CC0 dedication. By submitting a pull request, you agree to comply with this waiver of copyright interest.

Source-code contributions are also made available under the Apache Software License, Version 2.0 or later, as described in [`LICENSE`](LICENSE).
