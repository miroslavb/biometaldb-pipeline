#!/bin/bash
# Persistent monitor + auto-deploy for the biometaldb 3D tail-recovery on hive T07.
# Launched as a detached systemd-run unit on THIS box, so it is NOT turn-scoped:
# it survives user messages and session end (unlike a run_in_background watcher).
# Flow: poll T07 -> on DONE: pull out/full, regen paginated review UI, deploy to
# viewer/full, repoint homepage, restart molserver, Telegram-notify. Failures notify too.
set -uo pipefail

ROOT=/root/biometaldb-3d
VIEWER=/root/.hermes-agent2/biometaldb
LOG=$ROOT/recover_pipeline.log
H=root@100.121.152.8
SSH="ssh -o StrictHostKeyChecking=no -o ConnectTimeout=12"
TOKEN=$(grep -E '^TELEGRAM_BOT_TOKEN=' /root/.claude-tg-bridge/.env | head -1 | cut -d= -f2- | tr -d ' "')
CHAT=54300857

log(){ echo "[$(date -u +%F_%T)] $*" >>"$LOG"; }
tg(){ curl -sS --max-time 60 -F "chat_id=$CHAT" -F "text=$1" \
      "https://api.telegram.org/bot${TOKEN}/sendMessage" >/dev/null 2>&1 || true; }

log "===== monitor armed (pid $$) ====="
tg "🛰️ biometaldb recovery monitor армирован (persistent daemon на боксе, не turn-based). Авто пул+деплой по завершении T07."

# ---- 1. wait for DONE; heartbeat every 30 min; bail if service dies without DONE ----
last_hb=0
while true; do
  st=$(timeout 25 $SSH $H "if [ -f $ROOT/out/full/DONE ]; then echo DONE; \
       elif ! systemctl is-active --quiet bm-recover.service; then echo DIED; \
       else echo RUN:\$(wc -l < $ROOT/out/full/records.jsonl); fi" 2>/dev/null)
  case "$st" in
    DONE)  log "DONE sentinel detected"; break ;;
    DIED)  log "bm-recover died without DONE"
           tg "❌ biometaldb recovery: служба bm-recover остановилась БЕЗ DONE — нужна проверка на T07."
           exit 2 ;;
    RUN:*) now=$(date +%s)
           if [ $((now-last_hb)) -ge 1800 ]; then
             tg "⏳ recovery: ${st#RUN:}/9414 готово"; last_hb=$now; log "hb ${st#RUN:}/9414"
           fi ;;
    *)     : ;;  # transient ssh failure -> retry next loop
  esac
  sleep 60
done
tg "✅ T07 recovery DONE — тяну out/full и деплою full-review…"

# ---- 2. pull out/full (tar stream; -C avoids the cwd-wipe bug); 2 attempts ----
pull_ok=0
for attempt in 1 2; do
  log "pull attempt $attempt"
  if $SSH $H "tar czf - -C $ROOT out/full" | tar xzf - -C $ROOT 2>>"$LOG"; then
    pull_ok=1; break
  fi
  log "pull attempt $attempt failed"; sleep 10
done
if [ "$pull_ok" != 1 ] || [ ! -f $ROOT/out/full/manifest.json ]; then
  tg "❌ pull out/full FAILED (см. $LOG)"; log "pull FAILED / no manifest"; exit 3
fi
nrec=$(wc -l < $ROOT/out/full/records.jsonl 2>/dev/null || echo 0)
log "pulled: $nrec records, manifest present"

# ---- 3. generate UI, deploy, repoint homepage, restart ----
log "generating full review UI (review_full_ui.py)"
if ! python3 $ROOT/recon3d/review_full_ui.py --out $ROOT/out/full >>"$LOG" 2>&1; then
  tg "❌ review_full_ui.py упал (см. $LOG)"; exit 4
fi
mkdir -p $VIEWER/viewer/full
cp -a $ROOT/out/full/review.html $ROOT/out/full/index.json $ROOT/out/full/manifest.json $VIEWER/viewer/full/ 2>>"$LOG"
for d in struct img archive; do
  [ -d $ROOT/out/full/$d ] && cp -a $ROOT/out/full/$d $VIEWER/viewer/full/ 2>>"$LOG"
done
log "viewer/full populated"

# repoint homepage link ir100 -> full (backup first; only that one link mentions ir100/review.html)
cp -n $VIEWER/mol_server.py $VIEWER/mol_server.py.bak.pre-fullrepoint 2>/dev/null || true
sed -i 's#/viewer/ir100/review.html#/viewer/full/review.html#' $VIEWER/mol_server.py
if grep -q '/viewer/full/review.html' $VIEWER/mol_server.py; then
  log "homepage repointed -> /viewer/full/review.html"
else
  tg "⚠️ repoint sed не сработал — проверь mol_server.py"
fi
systemctl restart biometaldb-molserver && log "molserver restarted" || tg "⚠️ molserver restart failed"

# ---- summary + sentinel ----
summary=$(python3 -c "import json;m=json.load(open('$ROOT/out/full/manifest.json'));print(f\"{m.get('n_ok')}/{m.get('n')} built · {m.get('n_isomers_total')} structures · {m.get('n_valid_struct')} valid\")" 2>/dev/null)
echo "$summary" > $ROOT/PIPELINE_DONE
log "PIPELINE DONE: $summary"
tg "🧪 full-review задеплоен → https://mol.biometal.xyz/viewer/full/review.html
$summary
Остался step 4 (gbrain page + timeline) — делается в сессии."
