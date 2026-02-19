#!/bin/bash

###update allensdk download function to accommodate longer timeouts and download all data
### === Step: increase allensdk HTTP stream timeout ===

ALLENSDK_API_FILE=".venv/lib/python3.9/site-packages/allensdk/api/api.py"
if [[ ! -f "$ALLENSDK_API_FILE" ]]; then
    echo "❌ allensdk api.py not found at:"
    echo "   $ALLENSDK_API_FILE"
    exit 1
fi

# Backup once
if [[ ! -f "${ALLENSDK_API_FILE}.bak" ]]; then
    cp "$ALLENSDK_API_FILE" "${ALLENSDK_API_FILE}.bak"
fi

# ---- Patch 1: stream_zip_directory_over_http ----
if grep -q "timeout=(30, 1000)" "$ALLENSDK_API_FILE" | grep -q "stream_zip_directory_over_http"; then
    echo "✅ stream_zip_directory_over_http already patched — skipping"
else
    sed -i \
        's/def stream_zip_directory_over_http(url, directory, members=None, timeout=(9.05, 31.1)):/def stream_zip_directory_over_http(url, directory, members=None, timeout=(30, 1000)): ###increase to longer timeout in case of wifi interruption etc./' \
        "$ALLENSDK_API_FILE"
    echo "✅ Patched stream_zip_directory_over_http"
fi

# ---- Patch 2: stream_file_over_http ----
if grep -q "def stream_file_over_http(url, file_path, timeout=(30, 1000))" "$ALLENSDK_API_FILE"; then
    echo "✅ stream_file_over_http already patched — skipping"
else
    sed -i \
        's/def stream_file_over_http(url, file_path, timeout=(9.05, 31.1)):/def stream_file_over_http(url, file_path, timeout=(30, 1000)):/g' \
        "$ALLENSDK_API_FILE"
    echo "✅ Patched stream_file_over_http"
fi