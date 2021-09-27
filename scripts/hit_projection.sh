



curl localhost:18080/graphics/ | python -m json.tool

curl -d '{"graphics": [{"type": "circle", "color": {"r": 1, "g": 1, "b": 0, "a": 1}, "center": {"x": 20, "y": 30}, "r": 10}]}' -X PUT localhost:18080/graphics/ | python -m json.tool


curl -d '{"shutdown": true}' -X PUT localhost:18080/crow/ | python -m json.tool
