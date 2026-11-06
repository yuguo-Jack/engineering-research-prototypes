"""Normalize localized MSVC include traces so Ninja tracks headers reliably."""
import re
import subprocess
import sys

process=subprocess.Popen(sys.argv[1:],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
for line in iter(process.stdout.readline,b''):
    try:
        text=line.decode('utf8')
    except UnicodeDecodeError:
        text=line.decode('mbcs',errors='replace')
    text=re.sub(r'^\s*注意:\s*包含文件:', 'Note: including file:', text)
    sys.stdout.buffer.write(text.encode('utf8'))
    sys.stdout.buffer.flush()
raise SystemExit(process.wait())
