#!/bin/zsh
set -euo pipefail

ROOT="${0:A:h}"
VENV_PY="$ROOT/.venv/bin/python"
VENV_PYINSTALLER="$ROOT/.venv/bin/pyinstaller"
[[ -x "$VENV_PY" && -x "$VENV_PYINSTALLER" ]] || { echo "Run python3 -m venv .venv and .venv/bin/pip install -r requirements.txt first" >&2; exit 1; }
mkdir -p "$ROOT/pyinstaller-cache" "$ROOT/swift-cache" "$ROOT/dist"
export PYINSTALLER_CONFIG_DIR="$ROOT/pyinstaller-cache"
"$VENV_PYINSTALLER" --noconfirm --onedir --name PrimerMapEngine \
  --distpath "$ROOT/builddist" --workpath "$ROOT/build" --specpath "$ROOT" \
  --collect-all primer3 --add-data "$ROOT/reference_cache:reference_cache" \
  "$ROOT/src/primer_map.py" > "$ROOT/pyinstaller_build.log" 2>&1
CLANG_MODULE_CACHE_PATH="$ROOT/swift-cache" /usr/bin/swiftc \
  -module-cache-path "$ROOT/swift-cache" -target arm64-apple-macos12.0 \
  "$ROOT/src/PrimerMapApp.swift" -o "$ROOT/PrimerMapGUI" \
  > "$ROOT/swift_build.log" 2>&1

ICON_SOURCE="$ROOT/assets/PrimerMap-logo-en.png"
"$VENV_PY" -c 'from PIL import Image; import sys; Image.open(sys.argv[1]).convert("RGBA").resize((1024, 1024), Image.Resampling.LANCZOS).save(sys.argv[2], format="ICNS")' \
  "$ICON_SOURCE" "$ROOT/assets/PrimerMap-en.icns"

APP="$ROOT/dist/PrimerMap.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp "$ROOT/PrimerMapGUI" "$APP/Contents/MacOS/PrimerMap"
ditto "$ROOT/builddist/PrimerMapEngine" "$APP/Contents/Resources/engine"
cp "$ROOT/assets/PrimerMap-en.icns" "$APP/Contents/Resources/PrimerMap.icns"
cp "$ICON_SOURCE" "$APP/Contents/Resources/PrimerMap-logo-en.png"
cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
<key>CFBundleName</key><string>PrimerMap</string>
<key>CFBundleDisplayName</key><string>PrimerMap</string>
<key>CFBundleIdentifier</key><string>org.boyuan.primermap</string>
<key>CFBundleVersion</key><string>1.1.0</string>
<key>CFBundleShortVersionString</key><string>1.1.0</string>
<key>CFBundleDevelopmentRegion</key><string>en</string>
<key>CFBundleExecutable</key><string>PrimerMap</string>
<key>CFBundleIconFile</key><string>PrimerMap.icns</string>
<key>LSMinimumSystemVersion</key><string>12.0</string>
<key>NSHighResolutionCapable</key><true/>
<key>NSHumanReadableCopyright</key><string>Boyuan Wang</string>
</dict></plist>
PLIST
codesign --force --deep --sign - "$APP" > "$ROOT/codesign.log" 2>&1
rm -f "$ROOT/dist/PrimerMap-English-mac-arm64.zip"
ditto -c -k --sequesterRsrc --keepParent "$APP" "$ROOT/dist/PrimerMap-English-mac-arm64.zip"
echo "$APP"
echo "$ROOT/dist/PrimerMap-English-mac-arm64.zip"
