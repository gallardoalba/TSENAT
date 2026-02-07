# CircleCI setup for TSENAT

This document describes the minimal steps to enable CircleCI for the repository and add the `GITHUB_TOKEN` environment variable used by the `deploy_docs` job.

1. Enable the project in CircleCI
   - Sign in at https://circleci.com with your GitHub account.
   - Go to **Projects** → **Add Projects** and find `gallardoalba/TSENAT`.
   - Click **Set Up Project** and choose **Use existing config** (CircleCI will use `.circleci/config.yml`).

2. Add `GITHUB_TOKEN` to CircleCI
   - Create a GitHub personal access token (PAT): GitHub → Settings → Developer settings → Personal access tokens.
   - Recommended scopes:
     - For public repos: `public_repo` (or use a fine‑grained token limited to this repository with `Contents: Read & Write`).
     - For private repos: `repo` or a fine‑grained token with repository write access.
   - Copy the token value (you won't be able to view it again).
   - In CircleCI: open the project → **Project Settings** → **Environment Variables** → **Add Environment Variable**.
     - Name: `GITHUB_TOKEN`
     - Value: the PAT you copied

3. Verify builds and deploys
   - Trigger a test by pushing a commit to the `stable` branch (CircleCI is configured to run on `stable`):

```bash
git checkout -b ci-test
git commit --allow-empty -m "CI test"
git push origin ci-test:stable
```

   - Monitor the pipeline at https://app.circleci.com/pipelines/github/gallardoalba/TSENAT.
   - If `deploy_docs` runs, it will attempt to push the generated `docs/` to the `gh-pages` branch using `GITHUB_TOKEN`.

4. Security notes
   - Never commit tokens to the repo. Use CircleCI environment variables or contexts.
   - Prefer fine‑grained tokens with the minimum required permissions.
   - Revoke the token in GitHub if it is ever exposed.

5. Optional: use an Organization Context
   - If you manage multiple repos, create a CircleCI Context (Organization Settings → Contexts) and store `GITHUB_TOKEN` there, then reference the context from `config.yml`.

If you want, I can also add a short note to `.github/README.md` explaining that CircleCI is the canonical CI for `stable`.
