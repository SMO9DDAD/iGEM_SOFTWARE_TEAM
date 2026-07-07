import os, subprocess, threading, queue, json, time, uuid
from flask import Flask, render_template, request, jsonify, Response, stream_with_context

app = Flask(__name__)

_sessions = {}
_lock     = threading.Lock()


def make_sid():
    return str(uuid.uuid4())[:8]


def run_script(sid, script_path, cwd):
    q   = queue.Queue()
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    proc = subprocess.Popen(
        ["python3", "-u", script_path],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=0,
        cwd=cwd,
        env=env,
    )

    with _lock:
        _sessions[sid] = {"q": q, "proc": proc}

    def _reader():
        """Llegeix caràcter a caràcter i envia línies + prompts."""
        buf = ""
        try:
            while True:
                ch = proc.stdout.read(1)
                if not ch:
                    break
                if ch == "\n":
                    q.put({"type": "log", "text": buf})
                    buf = ""
                else:
                    buf += ch
                    # Envia el buffer acumulat com a "prompt" si sembla
                    # que el procés està esperant input (no arriba \n)
                    # Ho fem enviant l'estat parcial cada vegada que
                    # el buffer canvia i acaba amb ": " o "? "
                    if buf.rstrip().endswith((":", "?")):
                        q.put({"type": "prompt", "text": buf})
                        buf = ""
            # Qualsevol resta
            if buf:
                q.put({"type": "log", "text": buf})
            proc.wait()
            q.put({"type": "done", "code": proc.returncode})
        except Exception as e:
            q.put({"type": "error", "text": str(e)})
            q.put({"type": "done", "code": -1})
        finally:
            with _lock:
                _sessions.pop(sid, None)

    threading.Thread(target=_reader, daemon=True).start()


# ── RUTES ─────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/run", methods=["POST"])
def api_run():
    data        = request.json or {}
    module      = data.get("module", "main")
    project_dir = data.get("project_dir", "").strip()

    if not project_dir or not os.path.isdir(project_dir):
        return jsonify({"error": f"Directori no trobat: {project_dir}"}), 400

    script_map = {
        "main":     "main_nou.py",
        "consulta": "main_consulta.py",
        "main4":    "main4_nou.py",
    }
    script_name = script_map.get(module)
    if not script_name:
        return jsonify({"error": f"Mòdul desconegut: {module}"}), 400

    script_path = os.path.join(project_dir, script_name)
    if not os.path.isfile(script_path):
        return jsonify({"error": f"Script no trobat: {script_path}"}), 400

    sid = make_sid()
    run_script(sid, script_path, project_dir)
    return jsonify({"session_id": sid})


@app.route("/api/stream/<sid>")
def api_stream(sid):
    def generate():
        deadline = time.time() + 7200
        while time.time() < deadline:
            with _lock:
                sess = _sessions.get(sid)
            if sess is None:
                yield "data: " + json.dumps({"type": "done", "code": 0}) + "\n\n"
                return
            try:
                msg = sess["q"].get(timeout=1.0)
                yield "data: " + json.dumps(msg) + "\n\n"
                if msg.get("type") == "done":
                    return
            except queue.Empty:
                yield "data: " + json.dumps({"type": "ping"}) + "\n\n"

    return Response(
        stream_with_context(generate()),
        mimetype="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.route("/api/input", methods=["POST"])
def api_input():
    data = request.json or {}
    sid  = data.get("session_id", "")
    text = data.get("text", "")

    with _lock:
        sess = _sessions.get(sid)
    if not sess:
        return jsonify({"error": "Sessió no trobada"}), 404

    try:
        sess["proc"].stdin.write(text + "\n")
        sess["proc"].stdin.flush()
        return jsonify({"ok": True})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/kill/<sid>", methods=["POST"])
def api_kill(sid):
    with _lock:
        sess = _sessions.get(sid)
    if sess:
        try:
            sess["proc"].terminate()
        except Exception:
            pass
    return jsonify({"ok": True})


@app.route("/api/check_dir", methods=["POST"])
def api_check_dir():
    data = request.json or {}
    d    = data.get("path", "").strip()
    if not d or not os.path.isdir(d):
        return jsonify({"ok": False, "error": "Directori no trobat"})
    scripts = {s: os.path.isfile(os.path.join(d, s))
               for s in ["main_nou.py", "main_consulta.py", "main4_nou.py"]}
    return jsonify({"ok": True, "scripts": scripts})


if __name__ == "__main__":
    app.run(debug=True, port=5051, threaded=True)