from flask import Flask, request, Response, send_file
import subprocess
import os
import json
import sys
import psutil

app = Flask(__name__)

def check_system_resources():
    cpu_percent = psutil.cpu_percent(interval=1)
    memory_percent = psutil.virtual_memory().percent
    
    print(f"DEBUG: Current CPU usage: {cpu_percent}%")
    print(f"DEBUG: Current memory usage: {memory_percent}%")
    
    if cpu_percent > 75 or memory_percent > 75:
        return False, {
            'error': 'System is currently under heavy load. Please try again later.',
            'cpu_usage': cpu_percent,
            'memory_usage': memory_percent
        }
    return True, {
        'cpu_usage': cpu_percent,
        'memory_usage': memory_percent
    }

# Stream progress updates, including real-time script output
def progress_stream(input_file, output_file):
    try:
        # Step 1: Notify file upload success
        yield f"data: {json.dumps({'message': 'File successfully uploaded'})}\n\n"
        print("DEBUG: File successfully uploaded")
        sys.stdout.flush()

        # Step 2: Run the R script from the specified directory
        current_user = os.getenv("USER") or os.getenv("USERNAME")
        script_directory = f"/home/{current_user}/PanAzimuthWebAPI/"
        r_command = f"Rscript ANNotate_R.R {input_file}"
        print(f"DEBUG: Running command: {r_command} in directory: {script_directory}")
        sys.stdout.flush()

        # Use Popen to stream output line by line
        process = subprocess.Popen(
            r_command,
            shell=True,
            cwd=script_directory,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Stream the script's stdout
        for line in process.stdout:
            yield f"data: {json.dumps({'output': line.strip()})}\n\n"
            print(f"DEBUG: {line.strip()}")  # Debug log for stdout
            sys.stdout.flush()

        # Stream the script's stderr
        for line in process.stderr:
            yield f"data: {json.dumps({'error': line.strip()})}\n\n"
            print(f"DEBUG: STDERR: {line.strip()}")  # Debug log for stderr
            sys.stdout.flush()

        process.wait()  # Wait for process to complete
        if process.returncode != 0:
            yield f"data: {json.dumps({'error': f'R script failed with exit code {process.returncode}'})}\n\n"
            print(f"DEBUG: R script failed with exit code {process.returncode}")
            sys.stdout.flush()
            return

        # Notify script execution success
        yield f"data: {json.dumps({'message': 'R script finished successfully'})}\n\n"
        print("DEBUG: R script finished successfully")
        sys.stdout.flush()

        # Step 3: Check for output file
        if not os.path.exists(output_file):
            yield f"data: {json.dumps({'error': 'Output file not generated'})}\n\n"
            print("DEBUG: Output file not generated")
            sys.stdout.flush()
            return

        # Notify output file is ready
        yield f"data: {json.dumps({'message': 'Output file ready for download'})}\n\n"
        print("DEBUG: Output file ready for download")
        sys.stdout.flush()

    except Exception as e:
        yield f"data: {json.dumps({'error': str(e)})}\n\n"
        print(f"DEBUG: Exception occurred: {e}")
        sys.stdout.flush()

# Main route to handle RDS upload and processing
@app.route('/process_rds', methods=['POST'])
def process_rds():
    # Check system resources first
    resources_ok, resource_info = check_system_resources()
    if not resources_ok:
        return Response("System is currently under heavy load. Please try again in a few minutes.", status=503)

    if 'file' not in request.files:
        return Response("No file part in request", status=400)

    file = request.files['file']
    if file.filename == '':
        return Response("No selected file", status=400)

    # Save the uploaded RDS file
    input_file = f"/tmp/{file.filename}"
    output_file = input_file.replace(".rds", "_ANN.rds")
    print(f"DEBUG: Saving file to: {input_file}")
    file.save(input_file)

    try:
        # Stream progress updates
        return Response(progress_stream(input_file, output_file), content_type='text/event-stream')
    except Exception as e:
        print(f"DEBUG: Error during processing: {e}")
        return Response(json.dumps({'error': str(e)}), content_type='application/json', status=500)

# Route to download the output file
@app.route('/download_output', methods=['GET'])
def download_output():
    output_file = request.args.get('output_file')
    if not output_file or not os.path.exists(output_file):
        return Response("Output file not found", status=404)
    
    print(f"DEBUG: Returning output file: {output_file}")
    return send_file(output_file, as_attachment=True)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
