from django.db import models
import uuid
import os

def upload_to(instance, filename):
    return f"uploads/{instance.id}/{filename}"

class AnalysisJob(models.Model):
    JOB_TYPES = (
        ('BULK', 'Bulk RNA-seq'),
        ('SCRNA', 'Single-cell RNA-seq'),
    )
    STATUS_CHOICES = (
        ('PENDING', 'Pending'),
        ('PROCESSING', 'Processing'),
        ('COMPLETED', 'Completed'),
        ('FAILED', 'Failed'),
    )

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)
    job_type = models.CharField(max_length=10, choices=JOB_TYPES, default='BULK')
    data_file = models.FileField(upload_to=upload_to)
    is_normalized = models.BooleanField(default=False)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='PENDING')
    result_dir = models.CharField(max_length=255, blank=True, null=True)
    error_message = models.TextField(blank=True, null=True)

    def __str__(self):
        return f"{self.job_type} - {self.id} ({self.status})"
