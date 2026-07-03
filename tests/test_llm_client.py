from imgt_app import llm_client


def test_disabled_returns_none(monkeypatch):
    monkeypatch.setattr(llm_client, "_base_url", lambda: "")
    assert llm_client.llm_available() is False
    assert llm_client.llm_json("sys", "user") is None


def test_parses_openai_json(monkeypatch):
    monkeypatch.setattr(llm_client, "_base_url", lambda: "http://fake")

    class _Resp:
        def raise_for_status(self):
            pass

        def json(self):
            return {"choices": [{"message": {"content": '{"intent": "dossier"}'}}]}

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def post(self, *a, **k):
            return _Resp()

    monkeypatch.setattr(llm_client.httpx, "Client", _Client)
    assert llm_client.llm_json("sys", "user") == {"intent": "dossier"}


def test_failure_returns_none(monkeypatch):
    monkeypatch.setattr(llm_client, "_base_url", lambda: "http://fake")

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def post(self, *a, **k):
            raise RuntimeError("boom")

    monkeypatch.setattr(llm_client.httpx, "Client", _Client)
    assert llm_client.llm_json("sys", "user") is None
