#base on latest python 
FROM python:latest

#add program
COPY worker.py ./

#install dependencies
RUN python -m pip install pika numpy uproot awkward vector aiohttp requests

#command to run program
CMD [ "python", "-u", "./worker.py"] 
