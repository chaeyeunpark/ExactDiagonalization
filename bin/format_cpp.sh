#!/usr/bin/env bash
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # build_utils
PROJECT_SOURCE_DIR="$(dirname "${SCRIPT_DIR}")"
FILE_ARR=( $("${PROJECT_SOURCE_DIR}/bin/cpp_files.py" --include-examples --include-tests | jq -r '.[] | .name') )

if [[ "$1" == "--check" ]]; then
	echo "Check format"
	clang-format-13 --dry-run -Werror -i ${FILE_ARR[@]/%/}
	exit $?
else
	echo "Formatiing ${FILE_ARR[@]}"
	clang-format-13 -i ${FILE_ARR[@]/%/}
fi
