program rref
 integer::R
 integer::C
 R=length(M(:,1))
 C=length(M(1,:))
 do j=1,R
    for i=j+1:R
        M(j,:)=-(M(i,j)/M(j,j))*M(j,:)
        M(i,:)=M(i,:)+M(j,:)
    end
end
for j=R:-1:2
    for i=j-1:-1:1
        M(j,:)=-(M(i,j)/M(j,j))*M(j,:)
        M(i,:)=M(i,:)+M(j,:)
    end
end
for i=1:R
    M(i,:)=1/M(i,i)*M(i,:)
end
for i=1:R
    fprintf('a%d = %f\n',i,M(i,C))
end
end
