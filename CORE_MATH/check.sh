#!/bin/bash

MAKE=make

FUN="${!#}"
ARGS=("${@:1:$#-1}")

FILE="$(echo src/*/*/"$FUN".c)"
ORIG_DIR="$(dirname "$FILE")"

if ! [ -d "$ORIG_DIR" ]; then
    echo "Could not find $FUN"
    exit 1
fi

TMP_DIR="$(mktemp -d -t core-math.XXXXXX)"

trap 'rm -rf "$TMP_DIR"' EXIT

DIR="$TMP_DIR/toto/$(basename "$ORIG_DIR")"

mkdir "$TMP_DIR/toto"
cp -a "$ORIG_DIR" "$ORIG_DIR/../support" "$TMP_DIR/toto"
cp -a "$ORIG_DIR/../../generic" "$TMP_DIR"

if [ -n "${ARGS[0]}" ]; then
    KIND="${ARGS[0]}"
    ARGS=("${ARGS[@]:1}")
else
    SIZE=${FILE#src/binary}
    SIZE=${SIZE%%/*}
    case "$SIZE" in
        32)
            KIND=--exhaustive
            ;;
        *)
            KIND=--worst
    esac
fi

case "$KIND" in
    --exhaustive)
        "$MAKE"  -C "$DIR" clean
        "$MAKE"  -C "$DIR" check_exhaustive

        if [ "${#ARGS[@]}" -eq 0 ]; then
            MODES=("--rndn" "--rndz" "--rndu" "--rndd")
        else
            MODES=("${ARGS[@]}")
        fi
        for MODE in "${MODES[@]}"; do
            echo "Running exhaustive check in $MODE mode..."
            "$DIR/check_exhaustive" "$MODE"
        done
        ;;
    --worst)
        "$MAKE"  -C "$DIR" clean
        "$MAKE"  -C "$DIR" check_worst
        if [ "${#ARGS[@]}" -eq 0 ]; then
            MODES=("--rndn" "--rndz" "--rndu" "--rndd")
        else
            MODES=("${ARGS[@]}")
        fi
        for MODE in "${MODES[@]}"; do
            echo "Running worst cases check in $MODE mode..."
            "$DIR/check_worst" "$MODE" < "${FILE%.c}.wc"
        done
        ;;
    --special)
        "$MAKE"  -C "$DIR" clean
        "$MAKE"  -C "$DIR" check_special
        if [ "${#ARGS[@]}" -eq 0 ]; then
            MODES=("--rndn" "--rndz" "--rndu" "--rndd")
        else
            MODES=("${ARGS[@]}")
        fi
        for MODE in "${MODES[@]}"; do
            echo "Running special checks in $MODE mode..."
            "$DIR/check_special" "$MODE"
        done
        ;;
    *)
        echo "Unrecognized command"
        exit 1
esac
