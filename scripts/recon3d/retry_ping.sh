#!/bin/bash
# One-shot progress line for the biometaldb retry-tail run on hive T07.
H=root@100.121.152.8
SSH="ssh -o StrictHostKeyChecking=no -o ConnectTimeout=12"
timeout 25 $SSH $H "cd /root/biometaldb-3d; \
  state=RUN; [ -f out/full/RETRY_DONE ] && state=DONE; [ -f out/full/RETRY_FAIL ] && state=FAIL; \
  if [ \$state = RUN ] && ! pgrep -f recon3d/retry_tail.py >/dev/null; then state=NOPROC; fi; \
  python3 -c \"
import json,collections
c=collections.Counter(); rec=0; n=0
try:
  for line in open('out/full/retry.jsonl'):
    n+=1; r=json.loads(line); rb=(r.get('retry_meta') or {}).get('recovered_by')
    if rb: rec+=1
    else: c[r.get('fail_reason','?')]+=1
except Exception: pass
br=' '.join(f'{k}={v}' for k,v in sorted(c.items()))
print(f'\$state {n}/353 · recovered {rec} · still-failed {n-rec} · [{br}]')
\"" 2>/dev/null || echo "ssh-hiccup (will retry)"
