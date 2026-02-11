#!/usr/bin/env bash
set -euo pipefail

# Build pkgdown site and deploy the contents of docs/ to the gh-pages branch.
# Usage: bash tools/deploy_ghpages.sh

ROOT="$(git rev-parse --show-toplevel)"
echo "Root: $ROOT"

# Ensure docs exists and is non-empty
if [ ! -d "$ROOT/docs" ]; then
  echo "ERROR: docs/ directory not found in project root. Build the site first." >&2
  exit 1
fi
if [ -z "$(ls -A "$ROOT/docs" 2>/dev/null)" ]; then
  echo "ERROR: docs/ is empty. Build the site first." >&2
  exit 1
fi

echo "1) Generate coverage report..."
if [ -d "$ROOT/coverage" ]; then
  echo "Coverage report already exists, skipping generation."
else
  echo "Coverage report not found, generating..."
  Rscript "$ROOT/tools/run_coverage.R"
fi

echo "2) Copy coverage folder..."
cp -r "$ROOT/coverage" "$ROOT/docs/coverage" 2>/dev/null || true

# Ensure the PDF manual (if present) is copied into docs/articles so it is
# available on the published site at /articles/TSENAT.pdf
echo "3) Ensure PDF vignette is present in docs/articles..."
if [ -f "$ROOT/inst/doc/TSENAT.pdf" ]; then
  mkdir -p "$ROOT/docs/articles"
  echo "Copying PDF manual to docs/articles/TSENAT.pdf..."
  cp -f "$ROOT/inst/doc/TSENAT.pdf" "$ROOT/docs/articles/TSENAT.pdf"
elif [ -f "$ROOT/inst/ext/TSENAT.pdf" ]; then
  mkdir -p "$ROOT/docs/articles"
  echo "Copying PDF manual from inst/ext to docs/articles/TSENAT.pdf..."
  cp -f "$ROOT/inst/ext/TSENAT.pdf" "$ROOT/docs/articles/TSENAT.pdf"
else
  echo "No PDF manual found in inst/doc/ or inst/ext/; proceeding without PDF."
fi

echo "4) Preparing temporary clone..."
TMPDIR=$(mktemp -d)
repo_dir="$TMPDIR/repo"

# Ensure we have a remote named 'origin'
if ! git -C "$ROOT" remote get-url origin >/dev/null 2>&1; then
  echo "ERROR: no 'origin' remote configured for repository at $ROOT" >&2
  rm -rf "$TMPDIR"
  exit 1
fi

ORIGIN_URL=$(git -C "$ROOT" remote get-url origin)
echo "Cloning from: $ORIGIN_URL"
git clone --quiet "$ORIGIN_URL" "$repo_dir"
cd "$repo_dir"

# Ensure cleanup on exit
cleanup() {
  local rc=$?
  cd "$ROOT" >/dev/null 2>&1 || true
  rm -rf "$TMPDIR"
  return $rc
}
trap cleanup EXIT

echo "5) Checking gh-pages branch on origin..."
if git ls-remote --exit-code --heads origin gh-pages >/dev/null 2>&1; then
  echo "Found origin/gh-pages, fetching into local branch 'gh-pages'"
  git fetch origin gh-pages:gh-pages
  git checkout gh-pages
else
  echo "No gh-pages branch on origin, creating an orphan local 'gh-pages' branch"
  git checkout --orphan gh-pages
  git rm -rf . >/dev/null 2>&1 || true
fi

echo "6) Copying docs/ into temporary repo"
# Use rsync to copy and mirror docs/ content into the worktree root
rsync -a --delete --exclude=.git "$ROOT/docs/" "$repo_dir/"

git add -A

if [ -z "$(git status --porcelain)" ]; then
  echo "No changes to deploy. Exiting."
else
  git commit -m "Deploy pkgdown site: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
  echo "Pushing to origin/gh-pages (force-with-lease)..."
  # prefer --force-with-lease to reduce accidental clobbering
  git push origin gh-pages --force-with-lease
fi

echo "Deployment completed."
