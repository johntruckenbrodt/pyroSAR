

def reach(current, goal, nogos=None):
    nogos = [] if nogos is None else nogos
    counter = 0
    while current < goal:
        steps = [x+current for x in range(1, 7)]
        steps = [x for x in steps if x not in nogos and x < goal]
        if len(steps) > 0:
            current += max(steps)
            counter += 1
        else:
            return None
    return counter


def explore(steps, ladders, snakes):
    change = False
    for pos in range(len(steps) - 1):
        if steps[pos] is not None:
            targets = [[x, x] for x in range(pos + 1, 101) if x not in [y[0] for y in ladders + snakes]]
            targets += [x for x in ladders + snakes if x[0] > pos + 1]
            for dest in targets:
                moves = reach(pos + 1, dest[0], [x[0] for x in ladders + snakes])
                if moves is not None:
                    n = steps[pos] + moves
                    if steps[dest[1] - 1] is None or n < steps[dest[1] - 1]:
                        steps[dest[1] - 1] = n
                        change = True
    return steps, change


t = int(raw_input().strip())
for i in range(t):
    n = int(raw_input().strip())
    ladders = [[int(x) for x in raw_input().strip().split(' ')] for y in range(n)]
    m = int(raw_input().strip())
    snakes = [[int(x) for x in raw_input().strip().split(' ')] for y in range(m)]
    steps = [0] + [None]*99
    change = True
    while change:
        steps, change = explore(steps, ladders, snakes)
    print steps
    print steps[-1] if steps[-1] is not None else -1
