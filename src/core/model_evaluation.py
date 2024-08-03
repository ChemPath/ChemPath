from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sqlalchemy.orm import Session
from models import MLModel

def evaluate_model(model_id: int, X_test, y_test, session: Session):
    model = session.query(MLModel).get(model_id)
    y_pred = model.predict(X_test)
    
    metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'precision': precision_score(y_test, y_pred, average='weighted'),
        'recall': recall_score(y_test, y_pred, average='weighted'),
        'f1': f1_score(y_test, y_pred, average='weighted')
    }
    
    for key, value in metrics.items():
        setattr(model, f'metric_{key}', value)
    session.commit()
    
    return metrics
