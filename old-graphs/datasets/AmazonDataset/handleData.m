function [total_categories,D]=handleData(Dat,C,cat)
close all;
D=Dat;
total_categories=0;
%[total_categories, D]=createDataMatrix(C);

%figure(1);
%manipulateCategoriesOfData(Dat,cat);
manipulateProductSimilarity(D);
end
function [total_categories, D]=createDataMatrix(C)
%Data Matrix
%rows: ID numbers
%Columns : ASIN, Title, Group, Salesrank, #Similar, ProductList of Similar, #Categories, List of Categories, #reviews,#downloaded,#avg_rating
index=1;
C=C{1};
endi=size(C,1);
sizeD=549550;
D=cell(sizeD,20);
tStart=tic;
%mapObj=container.Map;
%contains all the different codes
%correspondig to different categories
total_categories=[];
for index=1:endi
    
    Sini=strsplit(C{index},' ');
    if(strcmp(Sini(1),'Id:'))
        id=str2double(Sini{2});
        
        if id==sizeD
            break
        end
        id=id+1;
        if(mod(id,1000)==0)
            id
        end
        index=index+1;
        Sini=strsplit(C{index},' ');
        if (strcmp(Sini(1),'ASIN:'))
            D{id,1}=(Sini{2});
            index=index+1;
            Sini=strsplit(C{index},' ');
            if(strcmp(Sini(1),'discontinued'))
                D{id,4}=0;
                D{id,5}=0;
                D{id,11}=0;
                D{id,10}=0;
                D{id,9}=0;
                continue;
            else(strcmp(Sini(1),'title:'));
                S=strjoin(Sini(2:end));
                D{id,2}=(S);
                index=index+1;
                Sini=strsplit(C{index},' ');
                if(strcmp(Sini(1),'group:'))
                    D{id,3}=(Sini{2});
                    index=index+1;
                    Sini=strsplit(C{index},' ');
                    
                    if(strcmp(Sini(1),'salesrank:'))
                        D{id,4}=str2double(Sini{2});
                        index=index+1;
                        Sini=strsplit(C{index},' ');
                        if(strcmp(Sini(1),'similar:'))
                            D{id,5}=str2double(Sini{2});
                            D{id,6}=strjoin(Sini(3:end));
                            index=index+1;
                            Sini=strsplit(C{index},' ');
                            if(strcmp(Sini(1),'categories:'))
                                D{id,7}=str2double(Sini{2});
                                %D{id,8}=['']
                                
                                %here I use matrix to store the codes
                                %for every category of each feature
                                %A=zeros((D{id,7}),10);
                                
                                %here I use vector to keep only the common
                                %categories between each feature
                                common_categories_vector=[];
                                for id1=1 : ((D{id,7}))
                                    if(strcmp(D{id,3},'DVD'))
                                        categ=findCategoriesForDVD(C{index+id1});
                                    else
                                        categ=findCategoriesForBooksAndMusic(C{index+id1});
                                    end
                                    total_categories=union(total_categories,categ);
                                    s=size(categ);
                                    %A(id1,1:s(2))=categ;
                                    common_categories_vector=union(common_categories_vector,categ);
                                end
                                %D{id,8}=A;
                                D{id,8}=common_categories_vector;
                                index=index+((D{id,7}))+1;
                                Sini=strsplit(C{index},' ');
                                if(strcmp(Sini(1),'reviews:'))
                                    
                                    D{id,9}=str2double(Sini{3});
                                    D{id,10}=str2double(Sini{5});
                                    D{id,11}=str2double(Sini{8});
                                end
                                
                            end
                        end
                        
                    end
                end
            end
        end
    end
end
tElapsed=toc(tStart)

%plotData(D);

end

function L=findCategoriesForBooksAndMusic(S)
Sini=strsplit(S,{'|','[',']'},'CollapseDelimiters',true);
L=str2double(Sini(3:2:length(Sini)));

before=size(L,2);
L=L(isfinite(L(1,: )));
after=size(L,2);
if(before-after~=0)
    %
    BandMdif=before-after;
end
end

function L=findCategoriesForDVD(S)
Sini=strsplit(S,{'[',']'},'CollapseDelimiters',true);
Leb=str2double(Sini(2:2:length(Sini)));

before=size(Leb,2);
L=Leb(isfinite(Leb(1,: )));
%
after=size(L,2);
if(before-after~=0)
    
    dvddif=before-after
end

end
function r=SimilarityFactor
end

function ind=find_similar(D,s)
s=strsplit(s,' ');
ind=find(ismember({D{:,1}}',s));

end
function manipulateProductSimilarity(D)
D=D(1:548552,:);

D=filterData(D);
figure(1);
hold on;
ind=0;
i=0;
while(i<=50)
ind=ind+1;
i=i+1;
list_id=find_similar(D,D{ind,6});
if (size(list_id,1)<=1)
 %keyboard;
    i=i-1;
 %keyboard
 continue;
end

    list_id=[2;list_id];
    
    x=D(list_id,:); 
    loc=(1:size(x,1));
    y=cell2mat(x(:,11));
    %stem(loc+15*i,mean(y)*ones(size(y)));
    %stem(loc+15*i,y);
    stem(loc+15*i,var(y)*ones(size(y)));
end
hold off;
end

function D=filterData(D)
%rated above 20 users maybe
D=D((([D{:,9}]>20)'),:);
end
function manipulateCategoriesOfData(D,cat)
% q=(strcmp(D(:,3),'Book'));
% w=([D{:,11}]~=0);
% I must further search between books of the same type to find the correlation
% Create function that can read the catecories field. ... . . :(
%
%must not include unranked data
%which codes to include in the simulation
%check if list empty

figure(1);
hold on;
%for k=1:10
category_codes=randsample(cat,1)';
category_cod=ones(size(D,1),size(category_codes,2));
for(i=1:size(category_codes,1))
    category_cod(:,i)=category_codes(i).*category_cod(:,i);
end
category_cod=num2cell(category_cod);
%keyboard

x=find_categories(category_cod,D);
if(x==0)
    %continue;
    return;
end
% ;
C=D((([D{:,11}]>1)'),:);
%C=D(((([D{:,9}]>20)')&([D{:,11}]>=1)')&(x),:);
%C=D(((([D{:,9}]>300)')&(strcmp(D(:,3),'Book'))&([D{:,11}]>=1)'),:);
datasetsize=size(C,1)
if((datasetsize==0)||(datasetsize<5))
    % continue;
    return;
end
%C=D((([D{:,11}]>=1)'),:);
%add also information about the number of reviews
%a big number of reviews with high rating better
y=[C{:,4}];
x=[C{:,11}];
%x1=0.009*[C{:,9}];
%x=x.*x1;
%ym=norm(y);
%xm=norm(x);

% figure(1);
% plot(y);
% title(sprintf('Average salesrank for this category= %g',mean(y)));

%x=x/xm;
%y=y/ym;


scatter(x',y');
X=[ones(size(x')),x'];
b=regress(y',X,0.7);
%xlim([0 1]);
%ylim([0 1]);
xaxis=linspace(0,5);
plot(xaxis,b(1)+b(2)*xaxis);
legend(sprintf('Products belonging to this category = %g',datasetsize));
xlabel('average rating') ;
ylabel('SalesRank') ;
%hold off;
%end


hold off;
end
function x= find_categories(category_cod,D)
AgtB = cellfun(@(x,y) ismember(x,y),category_cod, {D{:,8}}','UniformOutput',false );
x=cell2mat(AgtB);
y=ones(size(D,1),1);

for i=1:size(x,2)
    y=y&x(:,i);
end
x=y;
end