#base on latest python
FROM python:latest

#Add program

COPY controller.py ./

#install dependencies
RUN python -m pip install pika numpy awkward vector uproot requests aiohttp atlasopenmagic matplotlib


#Command to run program
CMD [ "python", "-u", "./controller.py"]
