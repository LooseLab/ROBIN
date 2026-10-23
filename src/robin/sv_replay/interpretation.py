"""Junction topology annotations; these do not establish complete SV events."""

def annotate(call):
    kind = call['kind']
    same = call['c1'] == call['c2']
    span = abs(call['p2'] - call['p1']) if same else None
    topology = kind.lower()
    if kind == 'BND':
        topology = 'interchromosomal' if not same else 'unresolved'
        if same and span:
            left, right = sorted([(call['p1'], call['o1']), (call['p2'], call['o2'])])
            signs = (left[1], right[1])
            if signs in [('+', '+'), ('-', '-')]:
                topology = 'inversion_like'
            elif signs == ('-', '+'):
                topology = 'duplication_like'
            elif signs == ('+', '-'):
                topology = 'deletion_like'
    return dict(topology=topology, interval_bp=span)
