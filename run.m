%%=============================uniform distribution for E70s====================================
Nclus = 10;
Ncopy = 1000;
SNR = 1/50;
numcoef = 1000;
shift_range = 0;
file_name = 'E70s.mat';

[ FBsPCA_data, z,zcost,info, Timing, imager, image0, ZA] = hMRA_uniform( Nclus, Ncopy, SNR, numcoef, shift_range, file_name);

% plot figures like Fig. 2
data = FBsPCA_data.data_raw;      % groundtruth images
datan = FBsPCA_data.data_sample;  % noisy observations
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.0 0.0], [0.0 0.0]);
figure;
j=[1:4];    % determine which image to plot
for i=1:4
    subplot(4,4,(i-1)*4+1);
    imagesc(data(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+2);
    imagesc(image0(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+3);
    imagesc(datan(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+4);
    imagesc(imager(:,:,j(i)));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

% relative estimation error
e1 = sqrt(sum(sum(sum((imager-image0(:,:,ZA)).^2)))/sum(sum(sum(image0.^2))));
disp(['relative error: ', num2str(e1)]);

% relative sPCA error
data = FBsPCA_data.data_raw; 
e2 = sqrt(sum(sum(sum((image0-data).^2)))/sum(sum(sum(data.^2))));
disp(['relative error of sPCA: ', num2str(e2)]);


%%=========================nonuniform distribution for E70s=========================================
Nclus = 10;
Ncopy = [500*ones(1,5), 1500*ones(1,5)];
SNR = 1/50;
numcoef = 1000;
file_name = 'E70s.mat';

[ FBsPCA_data, z, zcost, info, Timing, imager, image0, ZA] = hMRA_nonuniform( Nclus, Ncopy, SNR, numcoef, file_name )

% plot figures like Fig. 4
data = FBsPCA_data.data_raw;      % groundtruth images
datan = FBsPCA_data.data_sample;  % noisy observations
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.0 0.0], [0.0 0.0]);
figure;
j=[find(ZA==6), find(ZA==7), find(ZA==1), find(ZA==2)];    % determine which image to plot
for i=1:4
    subplot(4,4,(i-1)*4+1);
    imagesc(data(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+2);
    imagesc(image0(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+3);
    imagesc(datan(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+4);
    imagesc(imager(:,:,j(i)));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

% relative estimation error for each image
e = zeros(1,10);
for i=1:10
    e(i) = sqrt(sum(sum((imager(:,:,find(ZA==i))-image0(:,:,i)).^2))/sum(sum(image0(:,:,i).^2)));
end
disp('relative error for each image: ')
disp(e);

% total variation error of the estimated distribution
pi_est = zeros(10,1);
pi_est(ZA) = z.c;
pi_acc = Ncopy/10000';
e_tv = sum(abs(pi_est-pi_acc))/2;
disp(['total variation error: ', num2str(e_tv)]);


%%=============================uniform distribution for TrpV1====================================
Nclus = 10;
Ncopy = 1000;
SNR = 1/50;
numcoef = 1000;
shift_range = 0;
file_name = 'emd8117.mat';

[ FBsPCA_data, z, zcost, info, Timing, imager, image0, ZA] = hMRA_uniform( Nclus, Ncopy, SNR, numcoef, shift_range, file_name);

% plot figures like Fig. 7
data = FBsPCA_data.data_raw;      % groundtruth images
datan = FBsPCA_data.data_sample;  % noisy observations
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.0 0.0], [0.0 0.0]);
figure;
j=[1:4];    % determine which image to plot
for i=1:4
    subplot(4,4,(i-1)*4+1);
    imagesc(data(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+2);
    imagesc(image0(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+3);
    imagesc(datan(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+4);
    imagesc(imager(:,:,j(i)));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

% relative estimation error
e1 = sqrt(sum(sum(sum((imager-image0(:,:,ZA)).^2)))/sum(sum(sum(image0.^2))));
disp(['relative error: ', num2str(e1)]);


%%===================nonuniform distribution for yeast mitochindrial ribosome======================
Nclus = 10;
Ncopy = [500*ones(1,5), 1500*ones(1,5)];
SNR = 1/50;
numcoef = 1000;
file_name = 'emd3551.mat';

[ FBsPCA_data, z, zcost, info, Timing, imager, image0, ZA] = hMRA_nonuniform( Nclus, Ncopy, SNR, numcoef, file_name )

% plot figures like Fig. 8
data = FBsPCA_data.data_raw;      % groundtruth images
datan = FBsPCA_data.data_sample;  % noisy observations
subplot = @(m,n,p) subtightplot(m, n, p, [0.0 0.01], [0.0 0.0], [0.0 0.0]);
figure;
j=[find(ZA==6), find(ZA==7), find(ZA==1), find(ZA==2)];    % determine which image to plot
for i=1:4
    subplot(4,4,(i-1)*4+1);
    imagesc(data(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+2);
    imagesc(image0(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+3);
    imagesc(datan(:,:,ZA(j(i))));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    subplot(4,4,(i-1)*4+4);
    imagesc(imager(:,:,j(i)));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

% relative estimation error for each image
e = zeros(1,10);
for i=1:10
    e(i) = sqrt(sum(sum((imager(:,:,find(ZA==i))-image0(:,:,i)).^2))/sum(sum(image0(:,:,i).^2)));
end
disp('relative error for each image: ')
disp(e);

% total variation error of the estimated distribution
pi_est = zeros(10,1);
pi_est(ZA) = z.c;
pi_acc = Ncopy/10000';
e_tv = sum(abs(pi_est-pi_acc))/2;
disp(['total variation error: ', num2str(e_tv)]);


%%======================================error and SNR=========================================
Nclus = 10;
Ncopy = 1000;
numcoef = 1000;
shift_range = 0;
file_name = 'E70s.mat';

% e(i) saves relative error for SNR 1/(10*i)
RES={};
i=0;
for SNR=[1/10,1/50,1/100,1/150,1/200,1/300]
    i=i+1;
    [ FBsPCA_data, z, zcost, info, Timing, imager, image0, ZA] = hMRA_uniform( Nclus, Ncopy, SNR, numcoef, shift_range, file_name);
    Res.FBsPCA_data = FBsPCA_data;
    Res.z = z;
    Res.zcost = zcost;
    Res.info = info;
    Res.Timing = Timing;
    Res.imager = imager;
    Res.image0 = image0;
    Res.ZA = ZA;
    RES{i} = Res;
end

% plot Fig. 3 in the paper
subplot = @(m,n,p) subtightplot(m, n, p, [0.05 0.01], [0.05 0.05], [0.05 0.05]);
figure;
J=[10,50,100,150,200,300];
for i=1:6
    data=RES{i}.FBsPCA_data.data_raw;
    image0=RES{i}.image0;
    imager=RES{i}.imager;
    datan=RES{i}.FBsPCA_data.data_sample;
    ZA=RES{i}.ZA;
    e1=sqrt(sum(sum(sum((image0-data).^2)))/sum(sum(sum(data.^2))));
    e2=sqrt(sum(sum(sum((imager-image0(:,:,ZA)).^2)))/sum(sum(sum(image0.^2))));
    
    j=find(ZA==2);
    
    subplot(4,6,i);
    imagesc(data(:,:,2));
    colormap(gray);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    subplot(4,6,6+i);
    imagesc(datan(:,:,2));
    colormap(gray);
    title(['1/',num2str(J(i))]);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    subplot(4,6,12+i);
    imagesc(image0(:,:,2));
    colormap(gray);
    title([num2str(e1*100,'%0.2f'),'%']);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    subplot(4,6,18+i);
    imagesc(imager(:,:,j));
    colormap(gray);
    title([num2str(e2*100,'%0.2f'),'%']);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end


%%=================================error with different pi====================================
Nclus = 2;
SNR = 1/50;
numcoef = 1000;
file_name = 'E70s.mat';

e_tv = zeros(1,5);   % total variation error of pi 
e1 = zeros(1,5);     % relative error for the first class
e2 = zeros(1,5);     % relative error for the second class
for i=1:5
    Ncopy = [i*1e3, 1e4-i*1e3];
    [ FBsPCA_data, z, zcost, info, Timing, imager, image0, ZA] = hMRA_nonuniform( Nclus, Ncopy, SNR, numcoef, file_name )
    
    % compute tv error
    pi_est = zeros(2,1);
    pi_est(ZA) = z.c;
    pi_acc = Ncopy/1e4';
    e_tv(i) = sum(abs(pi_est-pi_acc))/2;

    % compute estimation errors of the two classes
    e1(i) = sqrt(sum(sum((imager(:,:,find(ZA==1))-image0(:,:,1)).^2))/sum(sum(image0(:,:,1).^2)));
    e2(i) = sqrt(sum(sum((imager(:,:,find(ZA==2))-image0(:,:,2)).^2))/sum(sum(image0(:,:,2).^2)));
end


%%===========================relative error with small shifts=================================
Nclus = 1;
Ncopy = 5000;
SNR = 1/50;
numcoef = 1000;
file_name = 'E70s.mat';

e = zeros(1,6);  % relative estimation errors under different shift amplitudes
for s=0:5
    shift_range = s;
    [ FBsPCA_data, z, zcost, info, Timing, imager, image0, ZA] = hMRA_uniform( Nclus, Ncopy, SNR, numcoef, shift_range, file_name);
    e(s+1) = sqrt(sum(sum(sum((imager-image0).^2)))/sum(sum(sum(image0.^2))));
end


%%=============================EM with uniform distribution===================================
Nclus = 10;
Ncopy = 1000;
SNR = 1/50;
numcoef = 1000;
shift_range = 0;
file_name = 'E70s.mat';

RES={};
for i=1:7
    Ndir = 2*2^i;
    [FBsPCA_data, Timing, res_em, res_op] = hMRA_EM_uniform(Nclus, Ncopy, Ndir, SNR, numcoef, shift_range, file_name);
    Res.FBsPCA_data = FBsPCA_data;
    Res.Timing = Timing;
    Res.res_em = res_em;
    Res.res_op = res_op;
    RES{i} = Res;
end

%compute estimation errors and plot Fig. 5 in the paper
e_em=zeros(7,10);
e_op=zeros(7,10);
t_em=zeros(7,1);
t_op=zeros(7,1);
t_inv=zeros(7,1);
ite_em=zeros(7,1);
for i=1:7
    image0=RES{i}.res_em.image0;
    imager=RES{i}.res_em.imager;
    ZA=RES{i}.res_em.ZA;
    t_em(i)=RES{i}.Timing.t_em;
    ite_em(i)=RES{i}.res_em.iter;
    for j=1:10
        e_em(i,j)=sqrt(sum(sum((imager(:,:,j)-image0(:,:,ZA(j))).^2))/sum(sum(image0(:,:,ZA(j)).^2)));
    end
    
    image0=RES{i}.res_op.image0;
    imager=RES{i}.res_op.imager;
    ZA=RES{i}.res_op.ZA;
    t_op(i)=RES{i}.Timing.t_opt;
    t_inv(i)=RES{i}.Timing.t_invariants;
    for j=1:10
        e_op(i,j)=sqrt(sum(sum((imager(:,:,j)-image0(:,:,ZA(j))).^2))/sum(sum(image0(:,:,ZA(j)).^2)));
    end
end

figure;
hold on;
for i=1:7
    plot(2*2^i*ones(1,10),e_em(i,1:10),'bo');
    plot(2*2^i*ones(1,10),e_op(i,1:10),'r^');
end
axis([2,512,0,0.5]);
set(gca, 'XScale', 'log');
legend('EM','MRA');
xticks([4,8,16,32,64,128,256]);
xlabel('Discretization of in-plane rotations');
ylabel('Relative error');

figure;
hold on;
plot([4,8,16,32,64,128,256]',t_em,'b-o');
plot([4,8,16,32,64,128,256]',t_op+t_inv,'r-^');
plot([4,8,16,32,64,128,256]',t_inv,'k-s');
for i=1:7
    text(2*2^i, t_em(i)/1.5, num2str(ite_em(i)));
end
axis([2,512,10,1e7]);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('EM','MRA','Invariants computing','Location','northwest');
xticks([4,8,16,32,64,128,256]);
xlabel('Discretization of in-plane rotations');
ylabel('CPU time')


%%============================EM with nonuniform distribution=================================
Nclus = 10;
SNR = 1/100;
numcoef = 1000;
shift_range = 0;
file_name = 'E70s.mat';

RES={};
for i=1:7
    Ndir = 2*2^i;
    [FBsPCA_data, Timing, res_em, res_op] = hMRA_EM_nonuniform(Nclus, Ncopy, Ndir, SNR, numcoef, shift_range, file_name);
    Res.FBsPCA_data = FBsPCA_data;
    Res.Timing = Timing;
    Res.res_em = res_em;
    Res.res_op = res_op;
    RES{i} = Res;
end

%compute estimation errors and plot Fig. 6 in the paper
e_em=zeros(7,10);
e_op=zeros(7,10);
t_em=zeros(7,1);
t_op=zeros(7,1);
t_inv=zeros(7,1);
ite_em=zeros(7,1);
for i=1:7
    image0=RES{i}.res_em.image0;
    imager=RES{i}.res_em.imager;
    ZA=RES{i}.res_em.ZA;
    t_em(i)=RES{i}.Timing.t_em;
    ite_em(i)=RES{i}.res_em.iter;
    for j=1:10
        e_em(i,j)=sqrt(sum(sum((imager(:,:,j)-image0(:,:,ZA(j))).^2))/sum(sum(image0(:,:,ZA(j)).^2)));
    end
    
    image0=RES{i}.res_op.image0;
    imager=RES{i}.res_op.imager;
    ZA=RES{i}.res_op.ZA;
    t_op(i)=RES{i}.Timing.t_opt;
    t_inv(i)=RES{i}.Timing.t_invariants;
    for j=1:10
        e_op(i,j)=sqrt(sum(sum((imager(:,:,j)-image0(:,:,ZA(j))).^2))/sum(sum(image0(:,:,ZA(j)).^2)));
    end
end

figure;
hold on;
for i=1:7
    plot(2*2^i*ones(1,10),e_em(i,1:10),'bo');
    plot(2*2^i*ones(1,10),e_op(i,1:10),'r^');
end
axis([2,512,0,0.5]);
set(gca, 'XScale', 'log');
legend('EM','MRA');
xticks([4,8,16,32,64,128,256]);
xlabel('Discretization of in-plane rotations');
ylabel('Relative error');

figure;
hold on;
plot([4,8,16,32,64,128,256]',t_em,'b-o');
plot([4,8,16,32,64,128,256]',t_op+t_inv,'r-^');
plot([4,8,16,32,64,128,256]',t_inv,'k-s');
for i=1:7
    text(2*2^i, t_em(i)/1.5, num2str(ite_em(i)));
end
axis([2,512,10,1e7]);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('EM','MRA','Invariants computing','Location','northwest');
xticks([4,8,16,32,64,128,256]);
xlabel('Discretization of in-plane rotations');
ylabel('CPU time')
































