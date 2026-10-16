#!/bin/bash
# Release dsXNAbuilder in two steps.
#
#   scripts/release.sh 0.2.0             # local only: checks, version, CHANGELOG, commit, tag,
#                                        # then the next -dev version
#   scripts/release.sh 0.2.0 --publish   # push main and the tag, and create the GitHub release
#
# The first step changes nothing outside this checkout: look at its two commits and the tag
# (git log -3; git show v0.2.0) before publishing. Zenodo archives every GitHub release (one DOI
# per version) once the repository is enabled in Zenodo, so enable it before publishing.
#
# To deploy the release to the web site, build from its tag (web/VERSION is then "0.2.0"):
#   git checkout v0.2.0 && scripts/build_web.sh && scripts/deploy_danslab.sh ../web-danslab
#   git checkout main
set -euo pipefail
cd "$(dirname "$0")/.."
VER=${1:?usage: scripts/release.sh X.Y.Z [--publish]}
MODE=${2:-}
TAG="v$VER"
PY=${PYTHON:-.venv/bin/python}
[ -x "$PY" ] || PY=python3
[[ "$VER" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "the version must be X.Y.Z, not $VER" >&2; exit 1; }
[ "$(git branch --show-current)" = main ] || { echo "release from main" >&2; exit 1; }

section() {  # body of the CHANGELOG section [$1]
  "$PY" - "$1" <<'PY'
import re, sys
text = open("CHANGELOG.md", encoding="utf-8").read()
m = re.search(r"^## \[" + re.escape(sys.argv[1]) + r"\][^\n]*\n(.*?)(?=^## \[|^\[[^\]]+\]: |\Z)",
              text, re.S | re.M)
sys.stdout.write(m.group(1).strip() + "\n" if m and m.group(1).strip() else "")
PY
}

# ---------------------------------------------------------------- step 2: publish
if [ "$MODE" = "--publish" ]; then
  git rev-parse -q --verify "refs/tags/$TAG" >/dev/null \
    || { echo "no tag $TAG: run scripts/release.sh $VER first" >&2; exit 1; }
  NOTES=$(mktemp)
  section "$VER" > "$NOTES"
  [ -s "$NOTES" ] || { echo "CHANGELOG.md has no section for $VER" >&2; exit 1; }
  git push origin main "$TAG"
  gh release create "$TAG" --title "dsXNAbuilder $VER" --notes-file "$NOTES"
  rm -f "$NOTES"
  echo "published $TAG"
  exit 0
fi
[ -z "$MODE" ] || { echo "unknown option: $MODE" >&2; exit 1; }

# ---------------------------------------------------------------- step 1: local release
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "commit your changes first" >&2; exit 1; }
git fetch -q origin
[ "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)" ] \
  || { echo "main and origin/main differ: pull or push first" >&2; exit 1; }
! git rev-parse -q --verify "refs/tags/$TAG" >/dev/null || { echo "tag $TAG exists" >&2; exit 1; }
[ -n "$(section Unreleased)" ] || { echo "CHANGELOG.md: write the changes under [Unreleased]" >&2; exit 1; }
"$PY" -m pytest -q

DATE=$(date +%F)
NEXT=$(echo "$VER" | awk -F. '{print $1 "." $2 + 1 ".0-dev"}')
"$PY" - "$VER" "$DATE" <<'PY'
import re, sys
ver, date = sys.argv[1:]


def sub(path, pattern, repl):
    text = open(path, encoding="utf-8").read()
    new, n = re.subn(pattern, repl, text, count=1, flags=re.M)
    if n != 1:
        raise SystemExit(f"{path}: no match for {pattern}")
    open(path, "w", encoding="utf-8").write(new)


sub("dsxnabuilder/__init__.py", r'^__version__ = "[^"]*"', f'__version__ = "{ver}"')
# the [Unreleased] changes become the release; an empty [Unreleased] starts again above them
sub("CHANGELOG.md", r"^## \[Unreleased\]\n", f"## [Unreleased]\n\n## [{ver}] - {date}\n")
with open("CHANGELOG.md", "a", encoding="utf-8") as fh:
    fh.write(f"\n[{ver}]: https://github.com/MMBAI-Lab/dsXNAbuilder/releases/tag/v{ver}\n")
sub("CITATION.cff", r'^version: .*$', f'version: "{ver}"')
text = re.sub(r"^date-released: .*\n", "", open("CITATION.cff", encoding="utf-8").read(), flags=re.M)
open("CITATION.cff", "w", encoding="utf-8").write(text.rstrip("\n") + f'\ndate-released: "{date}"\n')
PY
git commit -q -am "Release $VER"
# verbatim: the CHANGELOG headings start with "#", which git would strip as comments
{ echo "dsXNAbuilder $VER"; echo; section "$VER"; } | git tag -a --cleanup=verbatim "$TAG" -F -

"$PY" - "$NEXT" <<'PY'
import re, sys
path = "dsxnabuilder/__init__.py"
text = open(path, encoding="utf-8").read()
open(path, "w", encoding="utf-8").write(
    re.sub(r'^__version__ = "[^"]*"', f'__version__ = "{sys.argv[1]}"', text, count=1, flags=re.M))
PY
git commit -q -am "Start $NEXT"

git log --oneline -3
echo
echo "Tagged $TAG; main is now $NEXT. Nothing has left this checkout."
echo "Check:   git show $TAG"
echo "Deploy:  git checkout $TAG && scripts/build_web.sh && scripts/deploy_danslab.sh ../web-danslab; git checkout main"
echo "Publish: scripts/release.sh $VER --publish"
