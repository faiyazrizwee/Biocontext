from rest_framework import serializers
from .models import AnalysisJob

class AnalysisJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalysisJob
        fields = '__all__'
        read_only_fields = ('id', 'created_at', 'status', 'result_dir', 'error_message')
