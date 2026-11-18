"""Finite public fixtures; no periodic geometry implementation.

python3 -B periodic_cases.py --fixtures
python3 -B periodic_cases.py --check
python3 -B periodic_cases.py --verify
python3 -B periodic_cases.py --demo
"""
import json
import subprocess
import sys
import tempfile
from fractions import Fraction as F
from pathlib import Path

import coverage as api
from coverage.region_validation import validate_region
from coverage.rings import region_contains

P = F(10)


def box(a, b, c, d):
    return [(a, b), (c, b), (c, d), (a, d)]


def cases():
    return [
        ("cross", [("A", box(9, 0, 11, 1), [])], F(2), 1, [0]),
        ("band", [("A", box(0, 0, 6, 1), []),
                  ("B", box(5, 0, 11, 1), [])], F(10), 1, [-1, 1]),
        ("point", [("A", box(9, 0, 10, 1), []),
                   ("B", box(0, 1, 1, 2), [])], F(2), 2, [0, 0]),
        ("edge", [("A", box(9, 0, 10, 1), []),
                  ("B", box(0, 0, 1, 1), [])], F(2), 1, [0]),
        ("hole", [("A", box(7, 0, 15, 4), [box(9, 1, 11, 3)]),
                  ("B", box(4, 0, 8, 4), [])], F(36), 1, [-1, 0, 1]),
        ("duplicates", [("A", box(9, 0, 11, 1), []),
                        ("B", list(reversed(box(-1, 0, 1, 1))), [])],
         F(2), 1, [0]),
        ("concave", [("A", [(8, 0), (12, 0), (12, 1), (9, 1),
                            (9, 3), (12, 3), (12, 4), (8, 4)], [])],
         F(10), 1, [0]),
        ("slant", [("A", [(9, 0), (11, 0), (10, 2)], [])],
         F(2), 1, [0]),
    ]


def fixtures():
    result = []
    for name, shapes, area, count, windings in cases():
        images = []
        edges = 0
        for _, outer, holes in shapes:
            region = validate_region(outer, holes)
            flat = [p for ring in region.rings for p in ring]
            assert max(p[0] for p in flat) - min(p[0] for p in flat) < P
            edges += sum(map(len, region.rings))
            for k in range(-2, 3):
                images.append(tuple(tuple((x + k * P, y) for x, y in ring)
                                    for ring in region.rings))
        assert edges <= 80
        if name != "slant":
            # Unit-cell count is exact here: all these fixtures are
            # integer-coordinate, axis-aligned shapes.
            actual = sum(any(region_contains((F(2*x+1, 2), F(2*y+1, 2)), rings)
                             for rings in images)
                         for x in range(10) for y in range(-1, 5))
            assert actual == area, (name, actual, area)
        else:
            # Single non-self-overlapping triangle: base 2, height 2.
            assert area == F(2) * F(2) / 2
        result.append({"case": name, "area": str(area),
                       "components": count, "windings": windings})
    print("PERIODIC_FIXTURES " + json.dumps(result, sort_keys=True))


def builder():
    value = getattr(api, "build_periodic_report", None)
    if not callable(value):
        raise AssertionError("periodic report feature is not implemented")
    return value


def prepare(shapes, batch="day", shift=0):
    members = []
    for name, outer, holes in shapes:
        def moved(ring):
            return [(F(x) + shift, F(y)) for x, y in ring]
        revision = api.put_region(name, moved(outer), [moved(h) for h in holes])
        members.append((name, revision))
    api.seal_batch(batch, members)
    return batch


def assert_value_error(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def check_report(report, case):
    name, _, expected_area, count, windings = case
    assert report["format"] == 2 and report["period"] == P
    assert report["measurements"]["total_area"] == expected_area
    assert len(report["components"]) == count
    assert sum((i["area"] for i in report["measurements"]["areas"]), F()) == expected_area
    assert len(report["measurements"]["areas"]) == count
    assert len({c["component_id"] for c in report["components"]}) == count
    originals = {(r["region_id"], r["revision"]): r["rings"]
                 for r in report["inputs"]}
    lengths = []
    seen_windings = []
    boundary_ids = set()
    for component, area in zip(report["components"], report["measurements"]["areas"]):
        assert area["component_id"] == component["component_id"] and area["area"] > 0
        for boundary in component["boundaries"]:
            bid = boundary["boundary_id"]
            assert bid not in boundary_ids
            boundary_ids.add(bid)
            start = boundary["start"]
            assert all(type(v) is F for v in start) and 0 <= start[0] < P
            p = list(start)
            squares = []
            assert type(boundary["winding"]) is int
            seen_windings.append(boundary["winding"])
            assert boundary["segments"]
            for index, segment in enumerate(boundary["segments"]):
                delta = segment["delta"]
                assert len(delta) == 2 and all(type(v) is F for v in delta)
                assert delta != [F(0), F(0)]
                z = [p[i] + delta[i] for i in range(2)]
                squares.append(sum((v*v for v in delta), F()))
                assert segment["sources"]
                expected_sample = api.sample_boundary(report["report_id"], bid, index, F(1, 3))
                assert expected_sample["point"] == [(p[0]+delta[0]/3) % P, p[1]+delta[1]/3]
                assert len(expected_sample["sources"]) == len(segment["sources"])
                for source, sampled in zip(segment["sources"], expected_sample["sources"]):
                    ring = originals[(source["region_id"], source["revision"])][source["ring_index"]]
                    ei = source["edge_index"]
                    a, b = ring[ei], ring[(ei+1) % len(ring)]
                    endpoints = [[a[i] + source[t]*(b[i]-a[i]) for i in range(2)]
                                 for t in ("t0", "t1")]
                    offsets = [p[0]-endpoints[0][0], z[0]-endpoints[1][0]]
                    assert offsets[0] == offsets[1] and (offsets[0]/P).denominator == 1
                    assert p[1] == endpoints[0][1] and z[1] == endpoints[1][1]
                    parameter = source["t0"] + (source["t1"]-source["t0"])/3
                    assert sampled["parameter"] == parameter
                    assert all(sampled[k] == v for k, v in source.items())
                if name == "duplicates":
                    assert {s["region_id"] for s in segment["sources"]} == {"A", "B"}
                if name == "band":
                    assert delta[1] == 0
                p = z
            assert p == [start[0] + boundary["winding"]*P, start[1]]
            lengths.append({"boundary_id": bid, "length_squared": squares})
    assert sorted(seen_windings) == windings
    assert lengths == report["measurements"]["boundary"]
    assert api.open_package(api.export_report(report["report_id"])) == report


def disk_chain(case):
    build = builder()
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "periodic.sqlite"
        with api.open_archive(path) as archive:
            prepare(case[1])
            report = build("periodic", "day", P)
            check_report(report, case)
            assert report == build("periodic", "day", 10)
            assert_value_error(lambda: build("periodic", "day", F(11)))
            assert_value_error(lambda: api.build_report("periodic", "day"))
            api.build_report("flat", "day")
            assert_value_error(lambda: build("flat", "day", P))
            for invalid in (0, -1, True, 10.0):
                assert_value_error(lambda invalid=invalid: build("invalid", "day", invalid))
            saved = api.export_report("periodic")
            boundary = report["components"][0]["boundaries"][0]["boundary_id"]
            sample = api.sample_boundary("periodic", boundary, 0, F(1, 3))
            damaged = json.loads(saved)
            damaged["report"]["components"][0]["boundaries"][0]["winding"] += 1
            assert_value_error(lambda: api.open_package(json.dumps(damaged)))
            report["components"].clear()
            report["inputs"].clear()
            api.put_region(case[1][0][0], box(0, 0, 1, 1))
            assert api.export_report("periodic") == saved
            def fail(name, connection):
                if name == "report.after_payload":
                    connection.execute("INSERT INTO missing_periodic_table VALUES(1)")
            archive.failure_hook = fail
            import sqlite3
            try:
                build("retry", "day", P)
            except sqlite3.OperationalError:
                pass
            else:
                raise AssertionError("publication failure was swallowed")
            for table in ("reports", "receipts"):
                assert archive.connection.execute(
                    "SELECT count(*) FROM " + table + " WHERE report_id='retry'"
                ).fetchone()[0] == 0
            archive.failure_hook = None
            assert build("retry", "day", P)["measurements"]["total_area"] == case[2]
        with api.open_archive(path):
            assert api.export_report("periodic") == saved
            assert api.sample_boundary("periodic", boundary, 0, F(1, 3)) == sample
            check_report(api.read_report("periodic"), case)
    return {"case": case[0], "area": str(case[2]), "components": case[3],
            "windings": case[4], "sample": [str(v) for v in sample["point"]],
            "offline_equal": api.open_package(saved)["measurements"]["total_area"] == case[2],
            "reopen_equal": True}


def check():
    build = builder()
    for case in cases():
        for shift in (0, -20):
            with api.open_archive(":memory:"):
                prepare(case[1], shift=shift)
                report = build("r", "day", P)
                check_report(report, case)
                assert report == build("r", "day", P)
    for case in cases():
        if case[0] in ("band", "hole"):
            disk_chain(case)
    print("PERIODIC_CHECK_OK")


def main(argv):
    if argv == ["--fixtures"]:
        fixtures()
    elif argv == ["--check"]:
        check()
    elif argv == ["--verify"]:
        completed = subprocess.run(
            [sys.executable, "-B", "-m", "unittest", "discover",
             "-s", "tests", "-p", "test_*.py", "-v"], check=False)
        if completed.returncode:
            raise SystemExit(completed.returncode)
        check()
    elif argv == ["--demo"]:
        observed = [disk_chain(c) for c in cases() if c[0] in ("band", "hole")]
        print("POINTWISE_SCENARIO " + json.dumps(observed, sort_keys=True))
    else:
        raise SystemExit("use --fixtures, --check, --verify or --demo")


if __name__ == "__main__":
    main(sys.argv[1:])
