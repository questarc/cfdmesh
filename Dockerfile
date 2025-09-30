# Use FEniCS stable image as base
FROM quay.io/fenicsproject/stable:2019.1.0

# Install Python and pip
RUN apt-get update && apt-get install -y python3-pip gmsh

# Upgrade pip
RUN pip3 install --upgrade pip

# Install Python dependencies
RUN pip3 install streamlit==1.31.0 numpy==1.26.4 plotly==5.24.1 meshio==5.3.5 pyvista==0.44.0

# Copy your app code
COPY . /app
WORKDIR /app

# Expose Streamlit port
EXPOSE 8501

# Run Streamlit
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
