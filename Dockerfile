FROM python:3.11

WORKDIR /app
COPY . .
RUN pip install -U pip
RUN pip install -r requirements.txt

CMD ["python", "src/merge_files.py"]
