CurrentKeys=`aws iam list-access-keys --user-name $AWS_USERNAME | jq --raw-output '.AccessKeyMetadata[] | .AccessKeyId'`
aws iam create-access-key --user-name zhuy | 