# app/celery_tasks.py
from celery import Celery

celery = Celery(__name__, broker='redis://redis:6379/0')

def make_celery(app):
    celery.conf.update(app.config)
    return celery

@celery.task(bind=True)
def process_peptides(self, file_content):
    try:
        # Your processing logic here
        # Example: Generating a report
        report_content = "Processed peptides report based on file content"

        return report_content  # Return the content that will be rendered in the HTML
    except Exception as e:
        raise self.retry(exc=e)
