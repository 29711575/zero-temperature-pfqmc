#!/usr/bin/env python
import argparse
import json
import math


def reject_constant(value):
    raise ValueError("non-standard JSON constant: {0}".format(value))


def check_finite(value, path="root"):
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        raise ValueError("non-finite value at {0}".format(path))
    if isinstance(value, dict):
        for key, item in value.items():
            check_finite(item, "{0}.{1}".format(path, key))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            check_finite(item, "{0}[{1}]".format(path, index))


parser = argparse.ArgumentParser()
parser.add_argument("--zero", action="store_true")
parser.add_argument("files", nargs="+")
args = parser.parse_args()

for name in args.files:
    with open(name, "r") as stream:
        value = json.load(stream, parse_constant=reject_constant)
    check_finite(value)
    if args.zero:
        assert value["status"] == "complete"
        assert value["average_sign"] == 0
        assert value["sign_reweighted_observables_status"] == "unresolved_zero_average_sign"
        for key in ("S_pi", "S_pi_dq", "R_cdw"):
            assert value[key] is None
print("strict_json_pass files={0} zero={1}".format(len(args.files), int(args.zero)))
