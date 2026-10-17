#!/bin/bash
# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later
# Copy the built web app into a local checkout of MMBAI-Lab/web-danslab and
# point the dsxnabuilder card at it. Does NOT commit or push: a push to main
# of web-danslab republishes the whole lab site, so review and push by hand.
#
#   scripts/build_web.sh
#   scripts/deploy_danslab.sh /path/to/web-danslab
set -euo pipefail
cd "$(dirname "$0")/.."
SITE=${1:?usage: scripts/deploy_danslab.sh /path/to/web-danslab}
[ -f "$SITE/data/apps.en.json" ] || { echo "not a web-danslab checkout: $SITE" >&2; exit 1; }
[ -f web/dsxnabuilder.zip ] || { echo "run scripts/build_web.sh first" >&2; exit 1; }

mkdir -p "$SITE/public/dsxnabuilder"
rsync -a --delete --exclude '.version' --exclude '.ngl-version' web/ "$SITE/public/dsxnabuilder/"

# Stamp the deployed copy with this build's version (written by build_web.sh): new asset
# URLs for every deployment, and PAGE_VERSION for app.js to detect a stale cached page.
VER=$(cat web/VERSION)
D="$SITE/public/dsxnabuilder"
sed -i -E "s@(href=\"style\.css)\?v=[^\"]*\"@\1?v=$VER\"@g; s@(src=\"app\.js)\?v=[^\"]*\"@\1?v=$VER\"@g" "$D/index.html"
sed -i -E "s@(new URL\('worker\.js)\?v=[^']*'@\1?v=$VER'@g; s@^const PAGE_VERSION = '[^']*';@const PAGE_VERSION = '$VER';@" "$D/app.js"
grep -q "PAGE_VERSION = '$VER'" "$D/app.js" && grep -q "app.js?v=$VER" "$D/index.html" \
  || { echo "version stamp failed" >&2; exit 1; }
echo "stamped $VER"

python3 - "$SITE" <<'PY'
import json, sys
site = sys.argv[1]
for lang in ("en", "es"):
    path = f"{site}/data/apps.{lang}.json"
    with open(path, encoding="utf-8") as fh:
        apps = json.load(fh)
    for app in apps:
        if app.get("slug") == "dsxnabuilder":
            app["url"] = "https://danslab.xyz/dsxnabuilder/"
            app["embed"] = True
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(apps, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print(f"updated {path}")
PY
du -sh "$SITE/public/dsxnabuilder"
echo "Now: cd $SITE && npm run build (check) && git diff --stat, then commit/push when ready."
