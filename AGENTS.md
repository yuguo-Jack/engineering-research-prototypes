# Repository maintenance

- Each projects/<id> is independent: no runtime imports from sibling projects or author-specific paths.
- Read README.md, REPRODUCE.md and RESULTS.md before modifying an existing result.
- Keep archived measurements immutable. Record fresh validation separately and name its scope.
- Do not describe a supplied reference comparison as same-machine speedup.
- Keep dependencies, original-input checksums, commands and resource requirements explicit.
- Preserve source attribution and license notices; do not invent an open-source grant.
- Keep private preparation, credentials, datasets, binaries, caches and invalid candidates out of Git.
- Generate the overview using python tools/catalog.py; check with python tools/check_repo.py.
- Do not create a remote, publish pages or push unless the user explicitly requests it.
