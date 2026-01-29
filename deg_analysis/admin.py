from django.contrib import admin
from .models import AnalysisJob

@admin.register(AnalysisJob)
class AnalysisJobAdmin(admin.ModelAdmin):
    list_display = ('id', 'job_type', 'is_normalized', 'status', 'created_at')
    list_filter = ('status', 'job_type', 'is_normalized')
    readonly_fields = ('id', 'created_at', 'result_dir', 'error_message')
