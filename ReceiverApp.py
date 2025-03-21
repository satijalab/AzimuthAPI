from flask import Flask, request, Response, send_file
import subprocess
import os
import json
import sys
import psutil
from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
import requests
import threading

app = Flask(__name__)

# Global counter for concurrent requests
request_counter = 0
request_lock = threading.Lock()

def increment_counter():
    global request_counter
    with request_lock:
        request_counter += 1
        current_count = request_counter
    return current_count

def decrement_counter():
    global request_counter
    with request_lock:
        request_counter -= 1
        current_count = request_counter
    return current_count

def setup_logger(ip_address):
    """Setup a logger for the current request"""
    # Create logs directory if it doesn't exist
    if not os.path.exists('logs'):
        os.makedirs('logs')
    
    # Create a unique log filename with IP and timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f"logs/api_{ip_address}_{timestamp}.log"
    
    # Create a new logger
    logger = logging.getLogger(f'request_{ip_address}')
    logger.setLevel(logging.DEBUG)
    
    # Remove any existing handlers
    logger.handlers = []
    
    # Create file handler
    file_handler = RotatingFileHandler(log_filename, maxBytes=1024*1024, backupCount=5)
    file_handler.setLevel(logging.DEBUG)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Add handler to logger
    logger.addHandler(file_handler)
    
    return logger

def should_filter_line(line):
    # List of patterns to filter out
    filter_patterns = [
        "Traceback (most recent call last):",
        "File \"/home/ec2-user/miniconda3/envs/AzimuthNN_min/lib/python3.12/site-packages/tensorflow/python/pywrap_tensorflow.py\"",
        "import ssl",
        "import _ssl",
        "ssl.py",
        "File \"/home/ec2-user/R/x86_64-pc-linux-gnu-library/4.4/reticulate/python/rpytools/loader.py\"",
        "_find_and_load_hook",
        "return _run_hook",
        "module = hook()",
        "_find_and_load(name, import_)",
        "ImportError: /lib64/libcrypto.so.3: version `OPENSSL_3.3.0' not found",
        "Warning: Failed to load ssl module. Continuing without ssl support.",
        "^^^^^^"
    ]
    
    return any(pattern in line for pattern in filter_patterns)

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
def progress_stream_error(resource_info, logger):
    try:
        message = f'System is under heavy load (CPU: {resource_info["cpu_usage"]}%, Memory: {resource_info["memory_usage"]}%). Please try again later.'
        yield f"data: {json.dumps({'message': message})}\n\n"
        logger.error(message)

    except Exception as e:
        error_message = str(e)
        yield f"data: {json.dumps({'error': error_message})}\n\n"
        logger.error(f"Exception occurred: {error_message}")

# Stream progress updates, including real-time script output
def progress_stream(input_file, output_file, logger):
    try:
        # Step 1: Notify file upload success
        yield f"data: {json.dumps({'message': 'File successfully uploaded'})}\n\n"
        logger.info("File successfully uploaded")

        # Step 2: Run the R script from the specified directory
        current_user = os.getenv("USER") or os.getenv("USERNAME")
        script_directory = f"/home/{current_user}/PanAzimuthWebAPI/"
        # This will only work if the package is installed in this particular direcoty.
        r_command = f"Rscript ANNotate_R.R {input_file}"
        logger.info(f"Running command: {r_command} in directory: {script_directory}")

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
            if not should_filter_line(line.strip()):
                yield f"data: {json.dumps({'output': line.strip()})}\n\n"
                logger.info(line.strip())

        # Stream the script's stderr
        for line in process.stderr:
            if not should_filter_line(line.strip()):
                yield f"data: {json.dumps({'error': line.strip()})}\n\n"
                logger.error(f"STDERR: {line.strip()}")

        process.wait()  # Wait for process to complete
        if process.returncode != 0:
            error_message = f'R script failed with exit code {process.returncode}'
            yield f"data: {json.dumps({'error': error_message})}\n\n"
            logger.error(error_message)
            return

        # Notify script execution success
        yield f"data: {json.dumps({'message': 'R script finished successfully'})}\n\n"
        logger.info("R script finished successfully")

        # Step 3: Check for output file
        if not os.path.exists(output_file):
            error_message = 'Output file not generated'
            yield f"data: {json.dumps({'error': error_message})}\n\n"
            logger.error(error_message)
            return

        # Notify output file is ready
        yield f"data: {json.dumps({'message': 'Output file ready for download'})}\n\n"
        logger.info("Output file ready for download")

    except Exception as e:
        error_message = str(e)
        yield f"data: {json.dumps({'error': error_message})}\n\n"
        logger.error(f"Exception occurred: {error_message}")

def is_private_ip(ip_address):
    """Check if an IP address is private"""
    private_patterns = [
        '10.',      # 10.0.0.0/8
        '172.16.', '172.17.', '172.18.', '172.19.', '172.20.', '172.21.', '172.22.', '172.23.',
        '172.24.', '172.25.', '172.26.', '172.27.', '172.28.', '172.29.', '172.30.', '172.31.',  # 172.16.0.0/12
        '192.168.',  # 192.168.0.0/16
    ]
    return any(ip_address.startswith(pattern) for pattern in private_patterns)

def get_location_info(ip_address):
    """Get location information for an IP address"""
    if is_private_ip(ip_address):
        network_segment = ip_address.split('.')[0]
        if network_segment == '10':
            return {'country': 'Internal Network', 'city': 'Local', 'region': '10.x.x.x Range'}
        elif network_segment == '172':
            return {'country': 'Internal Network', 'city': 'Local', 'region': '172.16-31.x.x Range'}
        elif network_segment == '192':
            return {'country': 'Internal Network', 'city': 'Local', 'region': '192.168.x.x Range'}
    
    try:
        response = requests.get(f'http://ip-api.com/json/{ip_address}')
        if response.status_code == 200:
            data = response.json()
            if data['status'] == 'success':
                return {
                    'country': data.get('country', 'Unknown'),
                    'city': data.get('city', 'Unknown'),
                    'region': data.get('regionName', 'Unknown')
                }
    except Exception as e:
        print(f"Error getting location info: {str(e)}")
    return {'country': 'Unknown', 'city': 'Unknown', 'region': 'Unknown'}

# Main route to handle RDS upload and processing
@app.route('/process_rds', methods=['POST'])
def process_rds():
    # Increment request counter and get current count
    current_requests = increment_counter()
    
    try:
        # Get client IP address
        ip_address = request.remote_addr
        
        # Get location information
        location = get_location_info(ip_address)
        print(f"New API call from IP: {ip_address}")
        print(f"Location: {location['city']}, {location['region']}, {location['country']}")
        print(f"Current concurrent requests: {current_requests}")
        
        # Setup logger for this request
        logger = setup_logger(ip_address)
        logger.info(f"New request from IP: {ip_address} ({location['city']}, {location['region']}, {location['country']})")
        logger.info(f"Concurrent requests: {current_requests}")
        
        # Check system resources first
        resources_ok, resource_info = check_system_resources()
        if not resources_ok:
            logger.warning(f"System resources check failed: CPU {resource_info['cpu_usage']}%, Memory {resource_info['memory_usage']}%")
            decrement_counter()  # Decrement counter before error return
            return Response(progress_stream_error(resource_info, logger), content_type='text/event-stream')

        if 'file' not in request.files:
            logger.error("No file part in request")
            decrement_counter()  # Decrement counter before error return
            return Response("No file part in request", status=400)

        file = request.files['file']
        if file.filename == '':
            logger.error("No selected file")
            decrement_counter()  # Decrement counter before error return
            return Response("No selected file", status=400)

        # Save the uploaded RDS file
        input_file = f"/tmp/{file.filename}"
        output_file = input_file.replace(".rds", "_ANN.rds")
        logger.info(f"Saving file to: {input_file}")
        file.save(input_file)

        def generate():
            try:
                for chunk in progress_stream(input_file, output_file, logger):
                    yield chunk
            finally:
                remaining_requests = decrement_counter()
                print(f"Request completed. Remaining concurrent requests: {remaining_requests}")
                logger.info(f"Request completed. Remaining concurrent requests: {remaining_requests}")

        return Response(generate(), content_type='text/event-stream')
        
    except Exception as e:
        logger.error(f"Error during processing: {e}")
        decrement_counter()  # Decrement counter before error return
        return Response(json.dumps({'error': str(e)}), content_type='application/json', status=500)

# Route to download the output file
@app.route('/download_output', methods=['GET'])
def download_output():
    # Get client IP address and setup logger
    ip_address = request.remote_addr
    logger = setup_logger(ip_address)
    logger.info(f"Download request from IP: {ip_address}")
    
    output_file = request.args.get('output_file')
    if not output_file or not os.path.exists(output_file):
        logger.error(f"Output file not found: {output_file}")
        return Response("Output file not found", status=404)
    
    try:
        # Create a copy of the file in logs directory with IP and timestamp
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = f"logs/output_{ip_address}_{timestamp}.rds"
        
        # Copy the file
        import shutil
        shutil.copy2(output_file, log_filename)
        logger.info(f"Copied output file to logs: {log_filename}")
        
        # Return the original file to the user
        logger.info(f"Returning output file to user: {output_file}")
        return send_file(output_file, as_attachment=True)
    except Exception as e:
        error_message = f"Error copying file to logs: {str(e)}"
        logger.error(error_message)
        return Response(error_message, status=500)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
