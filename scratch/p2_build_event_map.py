"""Rule-based P2 freshwater event map (best effort; grid delegated 2026-11-24: "do your best with
what you've got"). Rules are fixed here before any event-level scoring:
 R1 exclude: Stage 0 overlap (FAZ84459), marine/mesocosm (CMO, Chesterfield), coastal Churchill
    (salinity unresolved), GENICE DNA (origin unresolved).
 R2 include only libraries with a collection date (microscape or SampleTown) AND coordinates.
    A library with a date but blank coordinates is rescued only if its SampleTown site+date matches
    exactly one coordinate pair among other libraries (logged as rescued).
 R3 event = campaign x station: single-linkage clusters of libraries whose coordinates are within
    2 km AND dates within 7 days. Station/site text is not used (inconsistent labels; WINRI clash).
 R2b coordinates: microscape lat/lon; if absent, SampleTown lat/lon. If both are present and
    differ by > 2 km, the library is excluded as a source conflict (Astra 251090518096134144).
 R4 everything else is excluded with its reason; nothing is imputed.
"""
import csv, math, datetime, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, 'results', 'p2', 'breadth')
R = list(csv.DictReader(open(os.path.join(D, 'event_adjudication_draft.tsv')), delimiter='\t'))
def fl(x):
    try: return float(x)
    except Exception: return None
EXCL = {'exclude_primary', 'exclude_freshwater'}
L = []
for r in R:
    L.append(dict(lib=r['library'], date=r['date'] or r['sampletown_date'], lat=fl(r['lat']), lon=fl(r['lon']),
                  slat=fl(r.get('sampletown_lat')), slon=fl(r.get('sampletown_lon')),
                  site=r['sampletown_site'], status=r['proposed_status'], reason=r['reason']))
def km(a, b):
    p1, p2 = math.radians(a['lat']), math.radians(b['lat']); dl = math.radians(b['lon'] - a['lon'])
    return 6371 * math.acos(min(1, math.sin(p1)*math.sin(p2) + math.cos(p1)*math.cos(p2)*math.cos(dl)))
for l in L:
    l['conflict'] = False
    if l['lat'] is not None and l['slat'] is not None:
        l['conflict'] = km(l, dict(lat=l['slat'], lon=l['slon'])) > 2.0
    elif l['lat'] is None and l['slat'] is not None:
        l['lat'], l['lon'] = l['slat'], l['slon']
for l in L:
    if l['status'] in EXCL: l['use'], l['why'] = False, 'R1 ' + l['reason']
    elif 'Churchill' in l['reason'] or 'GENICE' in l['reason']: l['use'], l['why'] = False, 'R1 ' + l['reason']
    elif not l['date']: l['use'], l['why'] = False, 'R4 no date'
    elif l['conflict']: l['use'], l['why'] = False, 'R2b microscape vs SampleTown coordinate conflict > 2 km'
    elif l['lat'] is None:
        cands = {(m['lat'], m['lon']) for m in L if m is not l and m['site'] and m['site'] == l['site']
                 and m['date'] == l['date'] and m['lat'] is not None}
        if l['site'] and len(cands) == 1:
            l['lat'], l['lon'] = cands.pop(); l['use'], l['why'] = True, 'R2 rescued coords via SampleTown site+date'
        else: l['use'], l['why'] = False, 'R4 no coordinates'
    else: l['use'], l['why'] = True, 'R2 dated with coordinates'
def km(a, b):
    p1, p2 = math.radians(a['lat']), math.radians(b['lat']); dl = math.radians(b['lon'] - a['lon'])
    return 6371 * math.acos(min(1, math.sin(p1)*math.sin(p2) + math.cos(p1)*math.cos(p2)*math.cos(dl)))
def dd(a, b):
    return abs((datetime.date.fromisoformat(a['date'][:10]) - datetime.date.fromisoformat(b['date'][:10])).days)
U = [l for l in L if l['use']]
parent = list(range(len(U)))
def find(i):
    while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i in range(len(U)):
    for j in range(i+1, len(U)):
        if km(U[i], U[j]) <= 2.0 and dd(U[i], U[j]) <= 7: parent[find(i)] = find(j)
groups = {}
for i, l in enumerate(U): groups.setdefault(find(i), []).append(l)
for k, (g, ls) in enumerate(sorted(groups.items(), key=lambda x: (min(l['date'] for l in x[1]), x[1][0]['lat']))):
    for l in ls: l['event'] = f"FW_E{k+1:03d}"
with open(os.path.join(D, 'event_map_best_effort.tsv'), 'w') as fh:
    fh.write('library\tevent\n')
    # runner requires all 232 libraries (run_p2_cohort.py:5868/3544); excluded rows get its skip value
    # CONTROL, with the reason kept in the audit file (Astra 251090680310841344)
    for l in L: fh.write(f"{l['lib']}\t{l.get('event', 'CONTROL') if l['use'] else 'CONTROL'}\n")
with open(os.path.join(D, 'event_map_best_effort_audit.tsv'), 'w') as fh:
    fh.write('library\tused\tevent\tdate\tlat\tlon\trule\n')
    for l in L: fh.write(f"{l['lib']}\t{l['use']}\t{l.get('event','')}\t{l['date']}\t{l['lat']}\t{l['lon']}\t{l['why']}\n")
import collections
print('libraries used', len(U), 'events', len(groups), 'sizes', sorted(collections.Counter(l['event'] for l in U).values(), reverse=True))
print(collections.Counter(l['why'].split(' ')[0] + ' ' + ' '.join(l['why'].split(' ')[1:3]) for l in L))
