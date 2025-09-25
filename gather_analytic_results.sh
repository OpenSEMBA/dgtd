#!/bin/bash

SRC_ROOT="./Exports"
DST_ROOT="./CollectedAnalytics"

mkdir -p "$DST_ROOT"

# Find all AnalyticRMS.txt files
find "$SRC_ROOT" -type f -name "AnalyticRMS.txt" | while read -r FILE; do
    # Extract the folder name that contains the resonance info
    PARENT_DIR=$(basename "$(dirname "$(dirname "$FILE")")")

    # Define new name
    NEW_NAME="AnalyticRMS_${PARENT_DIR}.txt"

    # Copy file into destination with the new name
    cp "$FILE" "$DST_ROOT/$NEW_NAME"
done
