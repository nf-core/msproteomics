#!/bin/bash
set -euo pipefail

# build_nfcore.sh — Build nf-fragpipe Docker container for nf-core pipelines
#
# Academic container (fcyucn/fragpipe:24.0 base + your licensed tools):
#   export MSFRAGGER_ZIP=/path/to/MSFragger-4.4.1.zip
#   export IONQUANT_ZIP=/path/to/IonQuant-1.11.20.zip
#   export DIATRACER_ZIP=/path/to/diatracer-2.2.1.zip
#   ./build_nfcore.sh academic
#
# Commercial container (FragPipePlus distribution):
#   export FRAGPIPEPLUS_ZIP=/path/to/FragPipePlus-24.0-linux.zip
#   ./build_nfcore.sh commercial
#
# Options:
#   TAG=myimage:1.0 ./build_nfcore.sh academic   # custom image tag

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-}"
TAG="${TAG:-fragpipe-nfcore:24.0}"

usage() {
    cat <<'EOF'
Usage: ./build_nfcore.sh <academic|commercial>

Academic container:
  Required env vars:
    MSFRAGGER_ZIP   Path to MSFragger zip (from https://msfragger-upgrader.nesvilab.org/upgrader/)
    IONQUANT_ZIP    Path to IonQuant zip (from https://msfragger-upgrader.nesvilab.org/ionquant/)
    DIATRACER_ZIP   Path to DiaTracer zip (from https://msfragger-upgrader.nesvilab.org/diatracer/)

Commercial container:
  Required env vars:
    FRAGPIPEPLUS_ZIP  Path to FragPipePlus zip (from https://www.fragmatics.com/)
    LICENSE_FILE      Path to commercial license.dat file

Optional:
  TAG=myimage:1.0   Override the Docker image tag (default: fragpipe-nfcore:24.0)
EOF
    exit 1
}

validate_file() {
    local var_name="$1"
    local file_path="$2"
    if [ ! -f "$file_path" ]; then
        echo "ERROR: $var_name file not found: $file_path" >&2
        exit 1
    fi
    # Verify it's a valid zip
    if ! unzip -t "$file_path" > /dev/null 2>&1; then
        echo "ERROR: $var_name is not a valid zip file: $file_path" >&2
        exit 1
    fi
    echo "  $var_name: $file_path (OK)"
}

build_academic() {
    : "${MSFRAGGER_ZIP:?Set MSFRAGGER_ZIP to path of MSFragger zip}"
    : "${IONQUANT_ZIP:?Set IONQUANT_ZIP to path of IonQuant zip}"
    : "${DIATRACER_ZIP:?Set DIATRACER_ZIP to path of DiaTracer zip}"

    echo "Building academic container..."
    echo "Validating inputs:"
    validate_file "MSFRAGGER_ZIP" "$MSFRAGGER_ZIP"
    validate_file "IONQUANT_ZIP" "$IONQUANT_ZIP"
    validate_file "DIATRACER_ZIP" "$DIATRACER_ZIP"

    BUILDDIR=$(mktemp -d)
    trap "rm -rf '$BUILDDIR'" EXIT

    echo "Extracting tool zips..."
    unzip -q "$MSFRAGGER_ZIP" -d "$BUILDDIR/fragpipe/"
    unzip -q "$IONQUANT_ZIP" -d "$BUILDDIR/fragpipe/"
    unzip -q "$DIATRACER_ZIP" -d "$BUILDDIR/fragpipe/"

    cp "$SCRIPT_DIR/Dockerfile.academic" "$BUILDDIR/Dockerfile"

    echo "Building Docker image: $TAG"
    docker build -t "$TAG" "$BUILDDIR"
}

build_commercial() {
    : "${FRAGPIPEPLUS_ZIP:?Set FRAGPIPEPLUS_ZIP to path of FragPipePlus zip}"
    : "${LICENSE_FILE:?Set LICENSE_FILE to path of commercial license.dat}"

    echo "Building commercial container..."
    echo "Validating inputs:"
    validate_file "FRAGPIPEPLUS_ZIP" "$FRAGPIPEPLUS_ZIP"
    if [ ! -f "$LICENSE_FILE" ]; then
        echo "ERROR: LICENSE_FILE not found: $LICENSE_FILE" >&2
        exit 1
    fi
    echo "  LICENSE_FILE: $LICENSE_FILE (OK)"

    BUILDDIR=$(mktemp -d)
    trap "rm -rf '$BUILDDIR'" EXIT

    echo "Extracting FragPipePlus distribution..."
    unzip -q "$FRAGPIPEPLUS_ZIP" -d "$BUILDDIR/fragpipe/"

    cp "$LICENSE_FILE" "$BUILDDIR/license.dat"
    cp "$SCRIPT_DIR/Dockerfile.commercial" "$BUILDDIR/Dockerfile"

    echo "Building Docker image: $TAG"
    docker build -t "$TAG" "$BUILDDIR"
}

case "$MODE" in
    academic)    build_academic ;;
    commercial)  build_commercial ;;
    *)           usage ;;
esac

echo ""
echo "SUCCESS: Built $TAG"
echo ""
echo "Verify with:"
echo "  docker run --rm $TAG msfragger --version"
echo "  docker run --rm $TAG ionquant --version"
echo "  docker run --rm $TAG diatracer --version"
echo "  docker run --rm $TAG philosopher version"
echo "  docker run --rm $TAG percolator --help 2>&1 | head -1"
