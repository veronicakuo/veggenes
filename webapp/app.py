"""veggenes web app — Flask backend serving the two-tab cancer-genomics explorer."""
import traceback
from flask import Flask, request, jsonify, send_from_directory

import veggenes_core as core

app = Flask(__name__, static_folder="static", static_url_path="")


@app.route("/")
def index():
    return send_from_directory("static", "index.html")


def _err(e):
    return jsonify({"error": str(e), "trace": traceback.format_exc()}), 500


@app.route("/api/studies")
def studies():
    try:
        return jsonify(core.find_studies(request.args.get("q", ""), limit=int(request.args.get("limit", 40))))
    except Exception as e:
        return _err(e)


@app.route("/api/cancer_types")
def cancer_types():
    try:
        return jsonify(core.list_cancer_types(request.args["study"]))
    except Exception as e:
        return _err(e)


@app.route("/api/frequencies")
def frequencies():
    try:
        study = request.args["study"]
        ct = request.args.get("cancer_type") or None
        top_n = int(request.args.get("top_n", 25))
        return jsonify(core.frequencies(study, ct, top_n))
    except Exception as e:
        return _err(e)


@app.route("/api/cooccurrence", methods=["POST"])
def cooccurrence():
    try:
        p = request.get_json(force=True)
        return jsonify(core.cooccurrence(
            p["study"], p.get("genes", []), p.get("cancer_type") or None,
            level=p.get("level", "sample"), alterations=p.get("alterations", "altered")))
    except Exception as e:
        return _err(e)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8050)
    ap.add_argument("--host", default="127.0.0.1")
    args = ap.parse_args()
    print(f"veggenes running at http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=False, threaded=True)
