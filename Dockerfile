# Use the official, slim Python 3.12 base image
FROM python:3.12-slim

# (Optional but good practice) Install system libraries for drawing
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcairo2 \
 && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy and install Python dependencies from requirements.txt
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of our application code
COPY . .

# Tell Render what port to listen on
ENV PORT=10000

# Use a shell to correctly start the server with the PORT variable
CMD ["sh", "-c", "uvicorn main:app --host 0.0.0.0 --port ${PORT:-10000}"]