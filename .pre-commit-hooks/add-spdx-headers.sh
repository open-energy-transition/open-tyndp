#!/bin/bash

# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

# Script to add Open-TYNDP SPDX header to files

HEADER_TEXT="SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>"

# Function to add Open-TYNDP header to a file
add_header() {
    local file="$1"
    
    # Skip if file already has Open-TYNDP header
    if grep -q "Contributors to Open-TYNDP" "$file"; then
        return 0
    fi
    
    # Only modify files that already have SPDX licensing information
    if ! grep -q "SPDX-License-Identifier" "$file"; then
        return 0
    fi
    
    # Create temporary file
    local temp_file
    temp_file=$(mktemp)
    
    # Detect comment style and insert accordingly
    
    if grep -q "^#.*SPDX" "$file"; then
        # Hash style: # comment (Python, YAML, TXT, SMK, shell, TOML, etc.)
        echo "# $HEADER_TEXT" > "$temp_file"
        cat "$file" >> "$temp_file"
    elif grep -q "SPDX-License-Identifier" "$file" && head -1 "$file" | grep -q "^\.\.\s*$"; then
        # RST style: .. comment block with indented content
        head -1 "$file" > "$temp_file"
        echo "  $HEADER_TEXT" >> "$temp_file"
        tail -n +2 "$file" >> "$temp_file"
    elif grep -q "^%.*SPDX" "$file"; then
        # LaTeX style: % comment
        echo "% $HEADER_TEXT" > "$temp_file"
        cat "$file" >> "$temp_file"
    elif grep -q "SPDX-License-Identifier" "$file" && head -1 "$file" | grep -q "^<!--"; then
        # HTML/Markdown style: <!-- comment -->
        head -1 "$file" > "$temp_file"
        echo "$HEADER_TEXT" >> "$temp_file"
        tail -n +2 "$file" >> "$temp_file"
    elif grep -q "SPDX-License-Identifier" "$file" && head -1 "$file" | grep -q "^@Comment{"; then
        # BibTeX style: @Comment{...}
        head -1 "$file" > "$temp_file"
        echo "$HEADER_TEXT" >> "$temp_file"
        tail -n +2 "$file" >> "$temp_file"
    elif grep -q "^REM.*SPDX" "$file"; then
        # Windows batch style: REM comment
        echo "REM $HEADER_TEXT" > "$temp_file"
        cat "$file" >> "$temp_file"
    elif [[ "$file" == *.ipynb ]] && grep -q "SPDX-License-Identifier" "$file"; then
        # Jupyter notebook style: JSON source arrays
        awk -v header="$HEADER_TEXT" '
        /SPDX-FileCopyrightText/ && !added {
            print "    \"" header "\\n\","
            print
            added = 1
            next
        }
        { print }
        ' "$file" > "$temp_file"
    else
        # Unknown comment style - skip to be safe
        rm "$temp_file"
        return 0
    fi
    
    # Replace original file
    mv "$temp_file" "$file"
    
    echo "Added Open-TYNDP SPDX header to $file"
}

# Process each file
for file in "$@"; do
    # Skip if file doesn't exist or is a directory
    if [[ ! -f "$file" ]]; then
        continue
    fi
    
    add_header "$file"
done