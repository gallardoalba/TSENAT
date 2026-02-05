Rscript scripts/build_pkgdown.R
set -e
REPO="https://github.com/gallardoalba/TSENAT.git"
TMPDIR=$(mktemp -d)
echo "tmpdir=$TMPDIR"
if git clone --branch gh-pages --single-branch "$REPO" "$TMPDIR"; then
  echo "Cloned existing gh-pages"
else
  echo "gh-pages missing; cloning repo and creating orphan"
  git clone "$REPO" "$TMPDIR"
  cd "$TMPDIR"
  git checkout --orphan gh-pages
fi
cd "$TMPDIR"
# remove all files except .git
find . -maxdepth 1 ! -name .git ! -name . -exec rm -rf {} +
# copy site
cp -R /home/nouser/galaxy/tools_source/TSENAT/docs/. . || true
# ensure there is something to commit
git add --all
if git diff --cached --quiet; then
  echo "No changes to push"
else
  git commit -m "Update pkgdown site from workspace docs" || true
  git push origin gh-pages
fi
printf "done\n"
echo continuing...; git --version; pwd; ls -la $TMPDIR; cd $TMPDIR || true; git status --porcelain --branch || true; printf 'done\n'
