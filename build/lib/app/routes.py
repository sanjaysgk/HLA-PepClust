from flask import Blueprint, request, jsonify, render_template
from .celery_tasks import process_peptides

bp = Blueprint("routes", __name__)


@bp.route("/upload")
def upload_form():
    return render_template("upload_form.html")


@bp.route("/process", methods=["POST"])
def process_file():
    file_content = request.files["file"].read().decode("utf-8")  # Read the file
    task = process_peptides.delay(file_content)  # Process the file with Celery
    return jsonify({"task_id": task.id}), 202


# New route to check task status and show result
@bp.route("/task_result/<task_id>")
def task_result(task_id):
    task = process_peptides.AsyncResult(
        task_id
    )  # Get the task using Celery's AsyncResult

    if task.state == "PENDING":
        # If the task is still pending, show a message
        content = "Task is still processing. Please check back later."
    elif task.state == "FAILURE":
        # If the task failed, show the exception traceback
        content = f"Task failed: {task.info}"
    elif task.state == "SUCCESS":
        # If the task is completed successfully, show the result
        content = task.result  # This should be the content you want to display (report)

    return render_template("report.html", content=content)
