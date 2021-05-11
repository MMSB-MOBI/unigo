from flask import Flask, request
import data_objects

def hello():
    data = request.get_json()
    schema = data_objects.CulledGoParametersSchema()
    parameters = schema.load(request)
    print(data)
    print(parameters)
    return 'hello'

app = Flask(__name__)
app.add_url_rule("/", 'hello', hello, methods=['POST'])

app.run(debug=True, port=3333)