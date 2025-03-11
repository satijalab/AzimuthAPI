from flask import Flask, request, send_file
import os

app = Flask(__name__)

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return 'No file part', 400

    file = request.files['file']
    if file.filename == '':
        return 'No selected file', 400

    # Save the uploaded file
    input_path = f'/tmp/{file.filename}'
    output_path = f'/tmp/uppercase_{file.filename}'
    file.save(input_path)

    # Convert to uppercase and save
    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            f_out.write(line.upper())

    # Return the file
    return send_file(output_path, as_attachment=True)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)