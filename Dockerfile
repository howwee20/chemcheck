# Use the official RDKit image as our base
FROM rdkit/rdkit:latest

# Set the working directory inside the container
WORKDIR /app

# Copy the requirements file and install our other dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of our application code into the container
COPY . .

# Tell Docker what command to run when the container starts
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "10000"]