%
%  FIGURES FOR THE PAPER ON MULTIKERNEL
%
%

classdef SemiParametricSimulations < simFunctionSet
	
	properties
		
	end
	
	methods
		%% Real data simulations 
		%  Data used: Swiss temperature
		% Goal: methods comparison NMSE
		function F = compute_fig_1001(obj,niter)
			
			%0. define parameters
			s_sigma=1.3;
			s_numberOfClusters=5;
			s_lambda=10^-5;
			s_monteCarloSimulations=niter;%100;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=1000;
			s_beta=0.02;
			s_alpha=0.005;
			s_niter=10;
			s_epsilon=0.2;
			s_functionTypes=5;
			v_sampleSetSize=(0.1:0.1:1);
			
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			% define graph
			[Ho,Mo,Alto,Hn,Mn,Altn] = readTemperatureDataset;
			%tic
			graphGenerator = GraphLearningSmoothSignalGraphGenerator('m_observed',Ho,'s_niter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta,'m_missingValuesIndicator',[]);
			%toc
			
			% tic
			% graphGenerator=SmoothSignalGraphGenerator('m_observed',Ho,'s_maxIter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta);
			% toc
			
			graph = graphGenerator.realization;
			%graph1 = graphGenerator1.realization;
			%L1=graph.getLaplacian
			%L2=graph1.getLaplacian
			v_sampleSetSize=round(v_sampleSetSize*graph.getNumberOfVertices);
			
			
			m_basis= SemiParametricSimulations.parametricPartForTempData(graph.getLaplacian,Alto,s_numberOfClusters);
			
			%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
			%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			%signal
			v_realSignal=Hn(:,3);
			functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignal,'s_normalize',1);
			% define bandlimited function estimator
			%m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth1);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth2);
			
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				
				
				s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
				%sample
				sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
				m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
				[m_samples,m_positions] = sampler.sample(m_graphFunction);
				%estimate
				m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
				m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
				m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
				
				% Performance assessment
				m_meanSquaredError(s_sampleSetIndex,1) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateBL1,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,2) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateBL2,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,3) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateNP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,4) = SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateSP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,5) = SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP,m_graphFunction);
				
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('real.mat');
			
			F = F_figure('X',v_sampleSetSize,'Y',m_meanSquaredError','xlab','Number of observed vertices (S)','ylab','NMSE','leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'});
		end
		%Goal: count only error on unseen data
		%     generalization capabibilities 
		function F = compute_fig_1011(obj,niter)
			
			%0. define parameters
			s_sigma=1.3;
			s_numberOfClusters=5;
			s_lambda=10^-5;
			s_monteCarloSimulations=niter;%100;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=1000;
			s_beta=0.02;
			s_alpha=0.005;
			s_niter=10;
			s_epsilon=0.2;
			s_functionTypes=5;
			v_sampleSetSize=(0.1:0.1:1);
		
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			% define graph
			[Ho,Mo,Alto,Hn,Mn,Altn] = readTemperatureDataset;
			%tic
			m_constraintLaplacian=zeros(size(Ho,1));
			m_constraintLaplacian(4:15,1)=1;
			graphGenerator = GraphLearningSmoothSignalGraphGenerator('m_observed',Ho,'s_niter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta,'m_constraintLaplacian',m_constraintLaplacian,'s_dont_estimate_the_signal',1);
			%toc
			
			% tic
			% graphGenerator=SmoothSignalGraphGenerator('m_observed',Ho,'s_maxIter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta);
			% toc
			
			graph = graphGenerator.realization;
			%graph1 = graphGenerator1.realization;
			%L1=graph.getLaplacian
			%L2=graph1.getLaplacian
			v_sampleSetSize=round(v_sampleSetSize*graph.getNumberOfVertices);
			
			
			m_basis= SemiParametricSimulations.parametricPartForTempData(graph.getLaplacian,Alto,s_numberOfClusters);
			
			%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
			%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			%signal
			v_realSignal=Hn(:,3);
			functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignal,'s_normalize',1);
			% define bandlimited function estimator
			%m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth1);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth2);
			
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				
				
				s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
				%sample
				sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
				m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
				[m_samples,m_positions] = sampler.sample(m_graphFunction);
				
				%estimate
				m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
				m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
				m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
				
				% Performance assessment
				m_indicator=SemiParametricSimulations.createIndicatorMatrix(m_graphFunction,m_positions);
				m_meanSquaredError(s_sampleSetIndex,1) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_indicator.*m_graphFunctionEstimateBL1,m_indicator.*m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,2) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_indicator.*m_graphFunctionEstimateBL2,m_indicator.*m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,3) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_indicator.*m_graphFunctionEstimateNP,m_indicator.*m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,4) = SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_indicator.*m_graphFunctionEstimateSP,m_indicator.*m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,5) = SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_indicator.*m_graphFunctionEpsilonInsesitiveEstimateSP,m_indicator.*m_graphFunction);
				
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('real.mat');
			
			F = F_figure('X',v_sampleSetSize,'Y',m_meanSquaredError','xlab','Number of observed vertices (S)','ylab','NMSE','leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'});
		end

		function F = compute_fig_1002(obj,niter)
			F = obj.load_F_structure(1001);
			F.ylimit=[0 1];
			F.xlimit=[10 89];
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			F.leg={strcat('Bandlimited  ',sprintf(' W=10')),strcat('Bandlimited',sprintf(' W=20')),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'};
			
			F.leg_pos = 'north';      % it can be 'northwest',
			%F.leg_pos_vec = [0.547 0.673 0.182 0.114];
		end
		function F = compute_fig_1012(obj,niter)
			F = obj.load_F_structure(1011);
			F.ylimit=[0 1];
			F.xlimit=[10 89];
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			F.leg={strcat('Bandlimited  ',sprintf(' W=10')),strcat('Bandlimited',sprintf(' W=20')),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'};
			
			F.leg_pos = 'north';      % it can be 'northwest',
			%F.leg_pos_vec = [0.547 0.673 0.182 0.114];
		end
		%% Synthetic data simulations 
		%  Data used: piece wise constant signals
		%  Goal: methods comparison NMSE
		function F = compute_fig_1003(obj,niter)
						
			%0. define parameters
			s_sigma=0.7;
			s_numberOfClusters=4;
			s_type=2;
			s_lambda=10^-8;
			s_monteCarloSimulations=niter;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=4;
			s_dataSetSize=100;
			s_functionTypes=5;
			s_epsilon=0.1;
			v_sampleSetSize=round((0.1:0.1:1)*s_dataSetSize);
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			% define graph function generator Parametric basis
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability',.6,'s_numberOfVertices',s_dataSetSize);
			graph = graphGenerator.realization;
			
			m_sparseBasis=graph.getClusters(s_numberOfClusters,s_type);
			m_basis=1.5*full(m_sparseBasis);
			%m_basis=m_basis*1;
			
			functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth1);
			functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			
			% define bandlimited function estimator
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth1);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth2);
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			%define semi-parametric epsinlon  insensitive estimator
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				
				
				s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
				%sample
				sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
				m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
				[m_samples,m_positions] = sampler.sample(m_graphFunction);
				%estimate
				m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
				m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
				m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
				% Performance assessment
				m_meanSquaredError(s_sampleSetIndex,1) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateBL1,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,2) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateBL2,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,3) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateNP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,4) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEstimateSP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,5) =SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP,m_graphFunction);
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('synthetic.mat');
			F = F_figure('X',v_sampleSetSize,'Y',m_meanSquaredError','xlab','Number of observed vertices (S)','ylab','NMSE','leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semiparametric (SL)','Semiparametric (\epsilon-IL)'});
			
			
		end
		
		function F = compute_fig_1004(obj,niter)
			F = obj.load_F_structure(1003);
			F.ylimit=[0 1];
			F.xlimit=[10 100];
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			F.leg={strcat('Bandlimited  ',sprintf(' W=10')),strcat('Bandlimited',sprintf(' W=20')),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'};
			
			%F.leg_pos = 'northeast';      % it can be 'northwest',
			F.leg_pos_vec = [0.647 0.683 0.182 0.114];
		end
		%% Real data simulations 
		%  Data used: Swiss temperature
		%  Goal: epsilon tuning
		function F = compute_fig_1005(obj,niter)
			%choose epsilon
			%0. define parameters
			s_sigma=1.3;
			s_numberOfClusters=5;
			s_lambda=10^-5;
			s_monteCarloSimulations=niter;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=1000;
			
			s_beta=0.02;
			s_alpha=0.005;
			s_niter=10;
			v_epsilonSet=[0,0.001,0.005,0.01,0.05,0.1,0.5,1,5];
			s_functionTypes=6;
			s_normalize=1;
			v_sampleSetSize=(0.5);
			
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			% define graph
			[Ho,Mo,Alto,Hn,Mn,Altn] = readTemperatureDataset;
			%tic
			graphGenerator = GraphLearningSmoothSignalGraphGenerator('m_observed',Ho,'s_niter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta,'m_missingValuesIndicator',[]);
			%toc
			
			% tic
			% graphGenerator=SmoothSignalGraphGenerator('m_observed',Ho,'s_maxIter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta);
			% toc
			
			graph = graphGenerator.realization;
			%graph1 = graphGenerator1.realization;
			%L1=graph.getLaplacian
			%L2=graph1.getLaplacian
			v_sampleSetSize=round(v_sampleSetSize*graph.getNumberOfVertices);
			
			
			m_basis= SemiParametricSimulations.parametricPartForTempData(graph.getLaplacian,Alto,s_numberOfClusters);
			
			%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
			%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			%signal
			v_realSignal=Hn(:,3);
			functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignal,'s_normalize',s_normalize);
			% define bandlimited function estimator
			%m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			%bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',m_laplacianEigenvectors(:,1:s_bandwidth1));
			%bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',m_laplacianEigenvectors(:,1:(s_bandwidth2)));
			
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			%nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			%semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				for s_epsilonSetIndex=1:size(v_epsilonSet,2)
					
					s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
					%sample
					sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
					m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
					[m_samples,m_positions] = sampler.sample(m_graphFunction);
					%estimate
					m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,v_epsilonSet(s_epsilonSetIndex));
					m_meanSquaredError(s_sampleSetIndex,s_epsilonSetIndex) = SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP,m_graphFunction);
				end
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('real.mat');
			
			F = F_figure('X',v_epsilonSet,'Y',m_meanSquaredError,'xlab','epsilon parameter','ylab','NMSE','leg',cellstr(num2str(v_sampleSetSize(:))),'logx',1);
		end
		
		%% Synthetic data simulations 
		%  Data used: piece wise constant signals
		% Goal: epsilon tuning
		function F = compute_fig_1007(obj,niter)
			%choose epsilon
			%0. define parameters
			s_sigma=1.3;
			s_numberOfClusters=5;
			s_lambda=10^-5;
			s_monteCarloSimulations=niter;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=6;
			s_type=2;
			s_beta=0.02;
			s_dataSetSize=100;
			s_alpha=0.005;
			s_niter=10;
			v_epsilonSet=[0,0.001,0.005,0.01,0.05,0.1,0.5,1,5];
			s_functionTypes=6;
			v_sampleSetSize=(0.5);
			
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability',.6,'s_numberOfVertices',s_dataSetSize);
			graph = graphGenerator.realization;
			
			m_sparseBasis=graph.getClusters(s_numberOfClusters,s_type);
			m_basis=1.5*full(m_sparseBasis);
			%m_basis=m_basis*1;
			
			functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth1);
			functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			
			v_sampleSetSize=round(v_sampleSetSize*graph.getNumberOfVertices);
			
			
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			%nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			%semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				for s_epsilonSetIndex=1:size(v_epsilonSet,2)
					
					s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
					%sample
					sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
					m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
					[m_samples,m_positions] = sampler.sample(m_graphFunction);
					%estimate
					m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,v_epsilonSet(s_epsilonSetIndex));
					m_meanSquaredError(s_sampleSetIndex,s_epsilonSetIndex) = SemiParametricSimulations.estimateNormalizedMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP,m_graphFunction);
				end
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('real.mat');
			
			F = F_figure('X',v_epsilonSet,'Y',m_meanSquaredError,'xlab','epsilon parameter','ylab','NMSE','leg',cellstr(num2str(v_sampleSetSize(:))),'logx',1);
		end
		%% Real data simulations 
		%  Data used: Swiss temperature
		%  Goal: Plot the gap for the epsilon insensitive function
		%  Data taken from CVX output just ploted here
		function F = compute_fig_1009(obj,niter)
			gap=[
				3.20E+04
				3.80E+02
				2.70E+01
				7.90E-01
				1.40E-02
				1.60E-04
				9.00E-06
				8.90E-07
				1.50E-08
				4.40E-10
				];
			F = F_figure('X',(1:size(gap,1)),'Y',gap','xlab','number of iterations','ylab','dual-primal gap','leg','Semi-parametric (\epsilon-IL)','logy',1);
			F.pos=[680 729 509 249];
		end
		
		
		
		
		%% Real data simulations
		%  Data used: Movielens Data
		%  Goal: compare the NMSE error for the known etries..
		%  IS NOT CORRECT VERSION DOES NOT FIND THE RIGHT NMSE
		% IMPLEMENT AGAIN HIDE ENTRIES
		function F=compute_fig_2001(obj,niter)
			%Movielens Simmulation
			%NMSE
			%0. define parameters
			s_sigma=0.12;
			s_lambda=10^-11;
			s_monteCarloSimulations=niter;
			s_bandwidth1=10;
			s_bandwidth2=5;
			s_fontSize=35;
			s_normalize=0;
			s_SNR=100;
			s_epsilon=0.5;
			s_functionTypes=5;
			s_type=2;
			v_sampleSetSize=(0.1:0.1:0.9);
			%specify number of items/signals
			s_totalItems=10;
			%# clusters for spectral
			s_numberOfClusters=5;
			
			m_meanSquaredErrorInner=zeros(size(v_sampleSetSize,2),s_functionTypes);
			m_meanSquaredErrorOuter=zeros(size(v_sampleSetSize,2),s_functionTypes);
			
			% define graph
			[m_clustUser,m_clustUserInfoAge,m_clustUserInfoSex,m_clustUserInfoOccup,m_train,m_test,c_userInfo,c_movieInfo]=prepareMLdat;
			m_adjacency= SemiParametricSimulations.generateAdjancency(m_train);
			graph = Graph('m_adjacency',m_adjacency);
			%define basis
			%m_basis=m_clustUserInfoOccup;
			m_basis=full(graph.getClusters(s_numberOfClusters,s_type));
			% define bandlimited function estimator
			m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth1);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth2);
			
			% define Kernel function
			%%Laplacian gives disconnected components must fix to proceeed maybe have
			%%different graphs.
			kernelgraphDiffusion = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix,'m_basis',m_basis);
			% Simulation
			
			for s_itemIndex=1:s_totalItems
				%Should sample only from the know values...
				
				v_realSignal=m_test(:,s_itemIndex);
				%indicator matrix so that only the known values are
				%considered for NMSE and signal estimation
				v_knownValuesInd=v_realSignal~=0;
				s_knownValues=nnz(v_realSignal);
				m_knownValuesInd=repmat(v_knownValuesInd, 1,s_monteCarloSimulations);
				v_sampleSetSizeForItem=round(v_sampleSetSize*s_knownValues);
				
				
				
				%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
				%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
				%signal
				functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignal,'s_normalize',s_normalize);
				
				
				
				
				for s_sampleSetIndex=1:size(v_sampleSetSizeForItem,2)
					
					s_numberOfSamples=v_sampleSetSizeForItem(s_sampleSetIndex);
					% % % %     if(s_numberOfSamples>s_monteCarloSimulations)
					% % % %         s_monteCarloSimulations=1;
					% % % %     end
					%sample
					sampler = PartiallyObservedGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR,'v_knownValuesInd',v_knownValuesInd);
					m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
					[m_samples,m_positions] = sampler.sample(m_graphFunction);
					%estimate
					m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
					m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
					m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
					m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
					m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
					
					% Performance assessment
					%compute the mse error only among the known entries
					m_meanSquaredErrorInner(s_sampleSetIndex,1) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateBL1.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,2) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateBL2.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,3) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateNP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,4) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateSP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,5) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
				end
				m_meanSquaredErrorOuter=m_meanSquaredErrorOuter+m_meanSquaredErrorInner;
			end
			m_meanSquaredErrorOuter=(1/s_totalItems)*m_meanSquaredErrorOuter;
			F = F_figure('X',100*v_sampleSetSize,'Y',m_meanSquaredErrorOuter','xlab','Percentage of observed vertices (S)','ylab','NMSE','tit',sprintf('#Items=%g',s_totalItems),'leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semiparametric (SL)','Semiparametric (\epsilon-IL)'});
			
			
		end
		function F = compute_fig_2002(obj,niter)
			F = obj.load_F_structure(2001);
			F.ylimit=[0 10];
			F.xlimit=[10 100];
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			
			%F.leg_pos = 'northeast';      % it can be 'northwest',
			F.leg_pos_vec = [0.647 0.683 0.182 0.114];
		end
		%% Real data simulations 
		%  Data used: Jensen joke Data
		%  Goal: compare the NMSE error for the known etries..
		%  WRONG VERSION DELETE GO TO 3003
		function F=compute_fig_3001(obj,niter)
			%Jensen Simmulation
			%NMSE
			%0. define parameters
			s_sigma=0.12;
			s_lambda=10^-11;
			s_monteCarloSimulations=niter;
			s_bandwidth1=10;
			s_bandwidth2=5;
			s_fontSize=35;
			s_normalize=0;
			s_SNR=100;
			s_epsilon=0.5;
			s_functionTypes=5;
			s_type=2;
			v_sampleSetSize=(0.8:0.1:0.9);
			%specify number of items/signals
			s_totalItems=3;
			%# clusters for spectral
			s_numberOfClusters=10;
			
			m_meanSquaredErrorInner=zeros(size(v_sampleSetSize,2),s_functionTypes);
			m_meanSquaredErrorOuter=zeros(size(v_sampleSetSize,2),s_functionTypes);
			
			% define graph
			m_reducedRatings=prepareJokerdat;
			%represent different maybe 0 rated jokes..
			m_reducedRatings((m_reducedRatings==99))=0;
			
			m_adjacency= SemiParametricSimulations.generateAdjancencyJensen(m_reducedRatings);
			graph = Graph('m_adjacency',m_adjacency);
			%define basis
			%m_basis=m_clustUserInfoOccup;
			m_basis=full(graph.getClusters(s_numberOfClusters,s_type));
			% define bandlimited function estimator
			m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth1);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth2);
			
			% define Kernel function
			%%Laplacian gives disconnected components must fix to proceeed maybe have
			%%different graphs.
			kernelgraphDiffusion = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix,'m_basis',m_basis);
			% Simulation
			
			for s_itemIndex=1:s_totalItems
				%Should sample only from the know values...
				
				v_realSignal=m_reducedRatings(:,s_itemIndex);
				%indicator matrix so that only the known values are
				%considered for NMSE and signal estimation
				v_knownValuesInd=v_realSignal~=0;
				s_knownValues=nnz(v_realSignal);
				m_knownValuesInd=repmat(v_knownValuesInd, 1,s_monteCarloSimulations);
				v_sampleSetSizeForItem=round(v_sampleSetSize*s_knownValues);
				
				
				
				%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
				%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
				%signal
				functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignal,'s_normalize',s_normalize);
				
				
				
				
				for s_sampleSetIndex=1:size(v_sampleSetSizeForItem,2)
					
					s_numberOfSamples=v_sampleSetSizeForItem(s_sampleSetIndex);
					% % % %     if(s_numberOfSamples>s_monteCarloSimulations)
					% % % %         s_monteCarloSimulations=1;
					% % % %     end
					%sample
					sampler = PartiallyObservedGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR,'v_knownValuesInd',v_knownValuesInd);
					m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
					[m_samples,m_positions] = sampler.sample(m_graphFunction);
					%estimate
					m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
					m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
					m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
					m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
					m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
					
					% Performance assessment
					%compute the mse error only among the known entries
					m_meanSquaredErrorInner(s_sampleSetIndex,1) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateBL1.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,2) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateBL2.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,3) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateNP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,4) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEstimateSP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					m_meanSquaredErrorInner(s_sampleSetIndex,5) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
				end
				m_meanSquaredErrorOuter=m_meanSquaredErrorOuter+m_meanSquaredErrorInner;
			end
			m_meanSquaredErrorOuter=(1/s_totalItems)*m_meanSquaredErrorOuter;
			F = F_figure('X',100*v_sampleSetSize,'Y',m_meanSquaredErrorOuter','xlab','Percentage of observed vertices (S)','ylab','NMSE','tit',sprintf('#users=%g',s_totalItems),'leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semiparametric (SL)','Semiparametric (\epsilon-IL)'});
			
			
		end
		
		function F = compute_fig_3002(obj,niter)
			F = obj.load_F_structure(3001);
			F.ylimit=[0 10];
			F.xlimit=[10 100];
			F.tit=sprintf('#users=%g',500);
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			
			%F.leg_pos = 'northeast';      % it can be 'northwest',
			F.leg_pos_vec = [0.647 0.683 0.182 0.114];
		end
		%% Real data simulations 
		%  Data used: Jensen joke Data
		%  Goal: compare the NMSE error excluding 2 jokes for each user
		function F=compute_fig_3003(obj,niter)
			%Jensen Simmulation
			%NMSE
			%excluding 2 jokes for each user
			%over these 2 I should do the monte carlos..
			%0. define parameters
			%TODO compute the RMSE online the test data..
			s_sigma=0.08;
			s_lambda=10^-11;
			s_monteCarloSimulations=niter;
			s_bandwidth1=2;
			s_bandwidth2=3;
			s_normalize=0;
			s_SNR=100;
			s_sigmaAdj=1;
			s_epsilon=0.5;
			s_functionTypes=5;
			s_type=3;
			%v_sampleSetSize=(0.8:0.1:0.9);
			
			%specify number of items/signals
			s_totalItems=100;
			%# clusters for spectral
			s_numberOfClusters=20;
			
			m_meanSquaredErrorInner=zeros(1,s_functionTypes);
			m_rootMeanSquaredErrorOuter=zeros(1,s_functionTypes);
			m_rootMeanSquaredErrorOuterOuter=zeros(1,s_functionTypes);
			% define graph
			m_reducedRatings=prepareJokerdat;
			%represent different maybe 0 rated jokes..
			m_reducedRatings((m_reducedRatings==99))=0;
			
			m_adjacency= SemiParametricSimulations.generateAdjancencyJensen(m_reducedRatings,s_sigmaAdj);
			graph = Graph('m_adjacency',m_adjacency);
			%define basis
			%m_basis=m_clustUserInfoOccup;
			m_basis=full(graph.getClusters(s_numberOfClusters,s_type));
			% define bandlimited function estimator
			m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth1);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',s_bandwidth2);
			
			% define Kernel function
			%%Laplacian gives disconnected components must fix to proceeed maybe have
			%%different graphs.
			kernelgraphDiffusion = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',kernelgraphDiffusion.generateKernelMatrix,'m_basis',m_basis);
			% Simulation
			%randomly exclude 2 measurments from each user
			m_reducedRatingsForEval=m_reducedRatings;
			m_unkwownInd=zeros(size(m_reducedRatings,1),2);
			for k=1:size(m_reducedRatings,1)
				v_ind=find(m_reducedRatings(k,:));
				m_unkwownInd(k,:)=datasample(v_ind,2,'Replace',false);
				m_reducedRatingsForEval(k,m_unkwownInd(k,:))=0;
			end
			s_numberOfTestValues=nnz(m_unkwownInd);
			for s_iterIndex=1:niter
				for s_itemIndex=1:s_totalItems
					%Should sample only from the know values...
					
					v_realSignalForEval=m_reducedRatingsForEval(:,s_itemIndex);
					v_realSignal=m_reducedRatings(:,s_itemIndex);
					%indicator matrix so that only the known values are
					%considered for NMSE and signal estimation
					v_knownValuesIndForEval=v_realSignalForEval~=0;
					v_knownValuesInd=v_realSignal~=0;
					
					s_knownValuesForEval=nnz(v_realSignalForEval);
					m_knownValuesInd=repmat(v_knownValuesInd, 1,s_monteCarloSimulations);
					%v_sampleSetSizeForItem=round(v_sampleSetSize*s_knownValues);
					
					
					
					%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
					%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
					%signal
					functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignalForEval,'s_normalize',s_normalize);
					
					
					
					
					s_sampleSetIndex=1;
					
					s_numberOfSamples=s_knownValuesForEval;
					% % % %     if(s_numberOfSamples>s_monteCarloSimulations)
					% % % %         s_monteCarloSimulations=1;
					% % % %     end
					%sample
					sampler = PartiallyObservedGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR,'v_knownValuesInd',v_knownValuesIndForEval);
					m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
					[m_samples,m_positions] = sampler.sample(m_graphFunction);
					%estimate
					m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
					m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
					m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
					m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
					%m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
					
					% Performance assessment
					%compute the mse error only among the known entries
					m_meanSquaredErrorInner(s_sampleSetIndex,1) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateBL1.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd);
					m_meanSquaredErrorInner(s_sampleSetIndex,2) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateBL2.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd);
					m_meanSquaredErrorInner(s_sampleSetIndex,3) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateNP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd);
					m_meanSquaredErrorInner(s_sampleSetIndex,4) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateSP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd);
					%m_meanSquaredErrorInner(s_sampleSetIndex,5) =SemiParametricSimulations.estimateRootMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP.*m_knownValuesInd,m_graphFunction.*m_knownValuesInd,s_knownValues);
					
					m_rootMeanSquaredErrorOuter=m_rootMeanSquaredErrorOuter+m_meanSquaredErrorInner;
				end
				%here evaluate the squareed
				m_rootMeanSquaredErrorOuter=sqrt((1/s_numberOfTestValues)*m_rootMeanSquaredErrorOuter);
				m_rootMeanSquaredErrorOuterOuter=m_rootMeanSquaredErrorOuterOuter+m_rootMeanSquaredErrorOuter;
			end
			m_rootMeanSquaredErrorOuterOuter=(1/niter)*m_rootMeanSquaredErrorOuterOuter;
			F = F_figure('plot_type_2D','bar','Y',m_rootMeanSquaredErrorOuterOuter,'xlab','Percentage of observed vertices (S)','ylab','NMSE','tit',sprintf('#jokes=%g',s_totalItems),'leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semiparametric (SL)','Semiparametric (\epsilon-IL)'});
			
			
		end
		function F = compute_fig_3004(obj,niter)
			F = obj.load_F_structure(3003);
			F.ylimit=[0 10];
			F.pos=[680 729 509 249];
			
			%F.leg_pos = 'northeast';      % it can be 'northwest',
			F.leg_pos_vec = [0.647 0.683 0.182 0.114];
		end
		%% Real data simulations 
		%  Data used: Temperature Dataset
		%  Goal: Convergence of objective function for the
		%  GraphLearningAlgorithm
		function F = compute_fig_1017(obj,niter)
			%convergence plot for estimate laplacian
			s_beta=0.002;
			s_alpha=0.02;
			s_niter=10;
			[Ho,Mo,Alto,Hn,Mn,Altn] = readTemperatureDataset;
			%tic
			graphGenerator = GraphLearningSmoothSignalGraphGenerator('m_observed',Ho,'s_niter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta,'m_missingValuesIndicator',[]);
			%toc
			
			% tic
			% graphGenerator=SmoothSignalGraphGenerator('m_observed',Ho,'s_maxIter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta);
			% toc
			
			[~,~,v_objective] = graphGenerator.realization;
			v_differenceObjective = -diff(v_objective);
			%s_omega=10^-6;
			%v_bound=s_omega*ones(size(v_differenceObjective));
			%Y=[v_differenceObjective;v_bound];
			F = F_figure('X',(2:size(v_objective,1)),'Y',v_differenceObjective,'xlab','number of iterations','ylab','objective function decrease','leg',{'GL-SigRep'});
		end
		function F = compute_fig_1018(obj,niter)
			F = obj.load_F_structure(1017);
			F.pos=[680 729 509 249];
            F.xlab='number of iterations';
            F.ylab='objective function decrease';
			F.logy=1;
			F.leg={'GL-SigRep'};
			%F.leg_pos = 'northeast';      % it can be 'northwest',
			F.leg_pos_vec = [0.647 0.683 0.182 0.114];
		end
		
		

		
	end
	
	
	
	
	methods(Static)
		%used for creating graph on recommendation systems maybe move....
		function m_adjacency=generateAdjancency(m_train)
			%% Create a graph of users based on their cosine similarity
			m_adjacency=zeros(size(m_train,1));
			for k=1:size(m_train,1)
				for l=1:k-1
					m_adjacency(k,l)=m_train(k,:)*(m_train(l,:))'/(norm(m_train(k,:))*norm(m_train(l,:)));
				end
			end
			m_adjacency(isnan(m_adjacency)) = 0 ;
			m_adjacency=m_adjacency+m_adjacency';
		end
		function m_adjacency=generateAdjancencyJensen(m_train,sigma)
			%% Create a graph of users based on their cosine similarity
			%Take into account negative values.. how to address..
			%as a first step take abs no abs just take the exp of the inner
			%product
			m_adjacency=zeros(size(m_train,1));
			for k=1:size(m_train,1)
				for l=1:k-1
					m_adjacency(k,l)=exp(sigma^2*m_train(k,:)*(m_train(l,:))'/(norm(m_train(k,:))*norm(m_train(l,:))));
				end
			end
			m_adjacency(isnan(m_adjacency)) = 0 ;
			m_adjacency=m_adjacency+m_adjacency';
		end
		function B=parametricPartForTempData(L,feat,n_clusters)
			
			
			C=kmeans(feat,n_clusters);
			% cluster via the altitude information..
			B=zeros(size(L,1),1);
			for i=1:n_clusters
				B(C==i,i)=1;
			end
			
		end
		
		function m_indicator=createIndicatorMatrix(m_graphFunction,m_samples)
			m_indicator=zeros(size(m_graphFunction));
			for i=size(m_graphFunction,2)
			m_indicator(m_samples(:,i),i)=1;
			end
			m_indicator=~m_indicator;
		end
		%estimates the normalized mean squared error
		function res = estimateNormalizedMeanSquaredError(m_est,m_observed)
			res=0;
			for i=1:size(m_est,2)
				res=res+norm(m_est(:,i)-m_observed(:,i))^2/norm(m_observed(:,i))^2;
			end
			res=(1/size(m_est,2))*res;
		end
		%estimates the root mean squared error over the known values
		function res = estimateMeanSquaredError(m_est,m_observed)
			res=0;
			for i=1:size(m_est,2)
				res=res+norm(m_est(:,i)-m_observed(:,i))^2;
			end
			res=(1/size(m_est,2))*res;
		end
		
		
		
		
	end
	
	
end





