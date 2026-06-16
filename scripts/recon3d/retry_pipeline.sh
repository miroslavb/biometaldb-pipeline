#!/bin/bash
# Persistent monitor + auto-deploy for the biometaldb FAILED-TAIL retry on hive T07.
# Launch detached on THIS box (systemd-run/setsid) so it survives turns & session end.
# Flow: poll T07 for RETRY_DONE -> pull out/full -> regen paginated review UI (now
# with failure-reason banners) -> deploy to viewer -> restart molserver -> TG-notify.
set -uo pipefail

ROOT=/root/biometaldb-3d
VIEWER=/root/.hermes-agent2/biometaldb
LOG=$ROOT/retry_pipeline.log
H=root@100.121.152.8
SSH="ssh -o StrictHostKeyChecking=no -o ConnectTimeout=12"
TOKEN=$(grep -E '^TELEGRAM_BOT_TOKEN=' /root/.claude-tg-bridge/.env | head -1 | cut -d= -f2- | tr -d ' "')
CHAT=54300857

log(){ echo "[$(date -u +%F_%T)] $*" >>"$LOG"; }
tg(){ curl -sS --max-time 60 -F "chat_id=$CHAT" -F "text=$1" \
      "https://api.telegram.org/bot${TOKEN}/sendMessage" >/dev/null 2>&1 || true; }

log "===== retry monitor armed (pid $$) ====="
tg "🛰️ biometaldb retry-tail монитор армирован (persistent daemon на боксе). Авто пул+деплой по завершении T07."

# ---- 1. wait for RETRY_DONE; heartbeat every 30 min; bail on RETRY_FAIL / dead proc ----
last_hb=0
while true; do
  st=$(timeout 25 $SSH $H "if [ -f $ROOT/out/full/RETRY_DONE ]; then echo DONE; \
       elif [ -f $ROOT/out/full/RETRY_FAIL ]; then echo FAIL; \
       elif ! pgrep -f retry_tail.py >/dev/null; then echo DEAD; \
       else echo RUN:\$(wc -l < $ROOT/out/full/retry.jsonl 2>/dev/null || echo 0); fi" 2>/dev/null)
  case "$st" in
    DONE) log "RETRY_DONE sentinel detected"; break ;;
    FAIL) log "RETRY_FAIL sentinel"; tg "❌ biometaldb retry: retry_tail вышел с ошибкой (см. out/retry_tail.log на T07)."; exit 2 ;;
    DEAD) log "retry_tail process gone without sentinel"
          tg "⚠️ biometaldb retry: процесс retry_tail исчез без DONE/FAIL — проверь T07."; exit 3 ;;
    RUN:*) now=$(date +%s)
           if [ $((now-last_hb)) -ge 1800 ]; then
             tg "⏳ retry: ${st#RUN:}/353 обработано"; last_hb=$now; log "hb ${st#RUN:}/353"
           fi ;;
    *) : ;;
  esac
  sleep 60
done
tg "✅ T07 retry DONE — тяну out/full и деплою обновлённый review…"

# ---- 2. pull out/full (tar stream; -C avoids cwd-wipe); 2 attempts ----
pull_ok=0
for attempt in 1 2; do
  log "pull attempt $attempt"
  if $SSH $H "tar czf - -C $ROOT out/full" | tar xzf - -C $ROOT 2>>"$LOG"; then
    pull_ok=1; break
  fi
  log "pull attempt $attempt failed"; sleep 10
done
if [ "$pull_ok" != 1 ] || [ ! -f $ROOT/out/full/manifest.json ]; then
  tg "❌ pull out/full FAILED (см. $LOG)"; log "pull FAILED / no manifest"; exit 4
fi

# ---- 3. regenerate UI (with failure-reason banners), deploy, restart ----
log "generating full review UI (review_full_ui.py)"
if ! python3 $ROOT/recon3d/review_full_ui.py --out $ROOT/out/full >>"$LOG" 2>&1; then
  tg "❌ review_full_ui.py упал (см. $LOG)"; exit 5
fi
mkdir -p $VIEWER/viewer/full
cp -a $ROOT/out/full/review.html $ROOT/out/full/index.json $ROOT/out/full/manifest.json $VIEWER/viewer/full/ 2>>"$LOG"
for d in struct img archive; do
  [ -d $ROOT/out/full/$d ] && cp -a $ROOT/out/full/$d $VIEWER/viewer/full/ 2>>"$LOG"
done
log "viewer/full populated"
systemctl restart biometaldb-molserver && log "molserver restarted" || tg "⚠️ molserver restart failed"

# ---- summary + sentinel ----
summary=$(python3 -c "import json;m=json.load(open('$ROOT/out/full/manifest.json'));print(f\"{m.get('n_ok')}/{m.get('n')} built (+{m.get('n_recovered',0)} recovered) · {m.get('n_still_failed',0)} still failed · reasons: {m.get('by_fail_reason')}\")" 2>/dev/null)
echo "$summary" > $ROOT/RETRY_DONE
log "RETRY PIPELINE DONE: $summary"
tg "🧪 retry-tail задеплоен → https://mol.biometal.xyz/viewer/full/review.html
$summary
Причины непостроенных видны в UI (фильтр why-failed + баннер на карточке)."
