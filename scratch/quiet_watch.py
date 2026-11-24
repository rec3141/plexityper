#!/usr/bin/env python3
"""Exit when #plexityper has had no post for QUIET seconds (default 600); print the last posts.
Read-only on the relay DB. Usage: quiet_watch.py [quiet_seconds]"""
import os, sqlite3, sys, time
QUIET = int(sys.argv[1]) if len(sys.argv) > 1 else 600
DB = os.path.expanduser('~/.agentworkforce/relay/local/relaycast.db')
Q = ("select m.id, a.name, m.created_at, substr(m.body,1,160) from messages m "
     "join channels ch on ch.id=m.channel_id join agents a on a.id=m.agent_id "
     "where ch.name='plexityper' order by m.created_at desc limit 6")
while True:
    c = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    rows = c.execute(Q).fetchall(); c.close()
    idle = time.time() - rows[0][2] if rows else QUIET
    if idle >= QUIET:
        print(f'#plexityper quiet for {int(idle)}s. Last posts:')
        for r in rows:
            print(r[0], r[1], time.strftime('%H:%M', time.localtime(r[2])), r[3].replace('\n', ' '))
        break
    time.sleep(max(30, QUIET - idle))
