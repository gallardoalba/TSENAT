# NA

CI alternatives and how to enable them

This repository includes example configuration files for several CI
providers so you can enable continuous integration without relying
solely on GitHub Actions minutes.

Files added: - `.circleci/config.yml` — CircleCI config (uses a rocker R
image). - `.travis.yml` — Travis CI config for R. - `.gitlab-ci.yml` —
GitLab CI config (uses a rocker R image).

How to enable each service from GitHub:

- CircleCI: Go to <https://circleci.com>, sign in with GitHub, enable
  the `gallardoalba/TSENAT` project. CircleCI will pick up
  `.circleci/config.yml` on push.
- Travis CI: Sign in at <https://travis-ci.com> with GitHub and enable
  this repository. Travis will run on pushes and pull requests using
  `.travis.yml`.
- GitLab CI: You can either import or mirror the GitHub repo into GitLab
  (<https://gitlab.com/new/import>) and enable CI; GitLab will run
  `.gitlab-ci.yml`. Alternatively, configure GitLab to mirror the GitHub
  repository so CI runs on GitLab when you push to GitHub.

Notes: - After enabling a service, you may need to grant access to the
repository and enable the CI service for the project. - These configs
are minimal; you may need to add OS packages or R packages specific to
your package’s dependencies. - If you want me to also add an
`appveyor.yml` or adapt configs for Windows builds, say so and I will
add them.
