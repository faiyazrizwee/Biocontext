from rest_framework import viewsets, status
from rest_framework.response import Response
from rest_framework.decorators import action
from .models import AnalysisJob
from .serializers import AnalysisJobSerializer
from .utils.analysis_engine import process_bulk_rna_seq, process_scrna_seq
import threading
import os

class AnalysisJobViewSet(viewsets.ModelViewSet):
    queryset = AnalysisJob.objects.all()
    serializer_class = AnalysisJobSerializer

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        job = serializer.save()
        
        # Start analysis in background
        thread = threading.Thread(target=self.run_analysis, args=(job.id,))
        thread.start()
        
        headers = self.get_success_headers(serializer.data)
        return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)

    def run_analysis(self, job_id):
        job = AnalysisJob.objects.get(id=job_id)
        job.status = 'PROCESSING'
        job.save()
        
        try:
            file_path = job.data_file.path
            
            if job.job_type == 'BULK':
                results = process_bulk_rna_seq(file_path, job.is_normalized)
            else:
                results = process_scrna_seq(file_path)
            
            # Save results to JSON file
            import json
            result_dir = f"media/results/{job.id}"
            os.makedirs(result_dir, exist_ok=True)
            result_file = f"{result_dir}/results.json"
            
            with open(result_file, 'w') as f:
                json.dump(results, f)
                
            job.result_dir = result_file
            job.status = 'COMPLETED'
            job.save()
            
        except Exception as e:
            job.status = 'FAILED'
            job.error_message = str(e)
            job.save()

    @action(detail=True, methods=['get'])
    def results(self, request, pk=None):
        job = self.get_object()
        if job.status == 'COMPLETED' and job.result_dir:
            import json
            try:
                with open(job.result_dir, 'r') as f:
                    data = json.load(f)
                return Response(data)
            except FileNotFoundError:
                return Response({'error': 'Result file not found.'}, status=status.HTTP_404_NOT_FOUND)
        elif job.status == 'FAILED':
             return Response({'error': job.error_message}, status=status.HTTP_400_BAD_REQUEST)
        else:
            return Response({'status': job.status}, status=status.HTTP_202_ACCEPTED)
