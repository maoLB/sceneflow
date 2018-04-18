function [Il_gray,Ir_gray,time_spent] = grayImages(Il,Ir)
% converts 4 views to grayscale images

tStart = tic;

n = numel(Il);

Il_gray = cell(1,n);
Ir_gray = cell(1,n);

for frame=1:n
    
    [~,~,c] = size(Il{frame});
    if c==3
        Il_gray{frame} = rgb2gray(Il{frame});
        Ir_gray{frame} = rgb2gray(Ir{frame});
    else
        Il_gray{frame} = Il{frame};
        Ir_gray{frame} = Ir{frame};
    end
    
end

time_spent = toc(tStart);
% fprintf('\nrgb2gray computed in\t\t%6.2f s\n',time_spent);

end