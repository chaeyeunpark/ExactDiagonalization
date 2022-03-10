#!/usr/bin/env bash
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # build_utils
PROJECT_SOURCE_DIR="$(dirname "${SCRIPT_DIR}")"
FILE_ARR=( $("${PROJECT_SOURCE_DIR}/bin/cpp_files.py" | jq -r '.[] | .name') )
echo "Formatiing ${FILE_ARR[@]}"
clang-format -i ${FILE_ARR[@]/%/}
