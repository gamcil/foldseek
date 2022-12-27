#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

DB="$1"
TMP_PATH="$2"

# shellcheck disable=SC2086
"$MMSEQS" mmcreateindex "${DB}" "${TMP_PATH}" ${CREATEINDEX_PAR} --index-subset 2 \
    || fail "createindex died"

# shellcheck disable=SC2086
"$MMSEQS" lndb "${DB}_h" "${DB}_ss_h" ${VERBOSITY_PAR} \
    || fail "lndb died"

# shellcheck disable=SC2086
"$MMSEQS" mmcreateindex "${DB}_ss" "${TMP_PATH}" ${CREATEINDEX_PAR} --index-subset 1 \
    || fail "createindex died"

# shellcheck disable=SC2086
"$MMSEQS" appenddbtoindex "${DB}_ca" "${DB}.idx" --id-list ${INDEX_DB_CA_KEY} ${VERBOSITY_PAR} \
    || fail "appenddbtoindex died"
